use bio_types::genome::{AbstractInterval, Interval, Locus};
use bio_types::strand::Strand;
use itertools::izip;
use itertools::Itertools;
use rust_htslib::faidx;

use crate::rada::modules::counting::{CountsBufferContent, LocusCounts, NucCounterContent};
use crate::rada::modules::dna::{Nucleotide, ReqNucleotide};
use crate::rada::modules::refnuc::RefNucPredictor;
use crate::rada::modules::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};

use super::{IntervalSummary, LocusSummary, MismatchSummary};

pub trait IntervalSummarizer {
    fn summarize(&self, name: &str, content: NucCounterContent) -> Vec<IntervalSummary>;
}

pub trait LocusSummarizer {
    fn summarize(&self, content: NucCounterContent) -> Vec<LocusSummary>;
}

pub trait Summarizer: IntervalSummarizer + LocusSummarizer {}

pub struct BaseSummarizer<RefNucPred: RefNucPredictor, StrandingPred> {
    refreader: faidx::Reader,
    refnucpred: RefNucPred,
    strandpred: StrandingPred,
}

impl<RefNucPred: RefNucPredictor, StrandingPred> BaseSummarizer<RefNucPred, StrandingPred> {
    pub(crate) fn new(
        refreader: faidx::Reader,
        refnucpred: RefNucPred,
        strandpred: StrandingPred,
    ) -> Self {
        BaseSummarizer {
            refreader,
            refnucpred,
            strandpred,
        }
    }

    fn reference(&self, interval: &Interval) -> Vec<Nucleotide> {
        let (start, end) = (
            interval.range().start as usize,
            interval.range().end as usize - 1,
        );
        let reference = self
            .refreader
            .fetch_seq(interval.contig(), start, end)
            .expect(&format!(
                "Failed to fetch sequence for region TODO:{}-{}",
                start, end
            ));
        assert_eq!(
            reference.len() as u64,
            interval.range().end - interval.range().start
        );

        reference.iter().map(|e| Nucleotide::from(*e)).collect()
    }

    fn sequence(
        &self,
        reference: &[Nucleotide],
        counts: &CountsBufferContent,
    ) -> Vec<ReqNucleotide> {
        let sumup = |a: &[LocusCounts], b: &[LocusCounts]| {
            a.iter().zip(b).map(|(a, b)| *a + *b).collect_vec()
        };
        let storage: Vec<LocusCounts>;
        let sequenced: &[LocusCounts] = match (counts.forward, counts.reverse, counts.unstranded) {
            (Some(f), None, None) => f,
            (None, Some(r), None) => r,
            (None, None, Some(u)) => u,
            (Some(f), Some(r), None) => {
                storage = sumup(f, r);
                &storage
            }
            (None, Some(r), Some(u)) => {
                storage = sumup(r, u);
                &storage
            }
            (Some(f), None, Some(u)) => {
                storage = sumup(f, u);
                &storage
            }
            (Some(f), Some(r), Some(u)) => {
                storage = izip!(f, r, u).map(|(a, b, c)| *a + *b + *c).collect_vec();
                &storage
            }
            (None, None, None) => {
                unreachable!()
            }
        };

        let variants = vec![];
        izip!(reference, sequenced)
            .map(|(refnuc, seqnuc)| self.refnucpred.predict(*refnuc, &variants, seqnuc))
            .collect()
    }

    fn mismatches(&self, sequence: &[ReqNucleotide], counts: &[LocusCounts]) -> MismatchSummary {
        let mut summary = MismatchSummary::zeros();
        for (seq, count) in sequence.iter().zip(counts) {
            match seq {
                ReqNucleotide::A => summary.A += *count,
                ReqNucleotide::C => summary.C += *count,
                ReqNucleotide::G => summary.G += *count,
                ReqNucleotide::T => summary.T += *count,
            }
        }
        summary
    }
}

impl<RefNucPred: RefNucPredictor, StrandingPred: IntervalStrandPredictor> IntervalSummarizer
    for BaseSummarizer<RefNucPred, StrandingPred>
{
    fn summarize(&self, name: &str, content: NucCounterContent) -> Vec<IntervalSummary> {
        // 1. get assembly data
        let reference = self.reference(&content.interval);

        // 2. predict reference nucleotides
        let counts = &content.counts;
        if counts.forward.is_none() && counts.reverse.is_none() && counts.unstranded.is_none() {
            // TODO: Log warning
            return vec![];
        }
        let sequence = self.sequence(&reference, counts);

        let mut result = Vec::with_capacity(3);
        if let Some(forward) = counts.forward {
            result.push(IntervalSummary {
                interval: content.interval.clone(),
                strand: Strand::Forward,
                name: name.to_owned(),
                mismatches: self.mismatches(&sequence, forward),
            })
        }

        if let Some(reverse) = counts.reverse {
            result.push(IntervalSummary {
                interval: content.interval.clone(),
                strand: Strand::Reverse,
                name: name.to_owned(),
                mismatches: self.mismatches(&sequence, reverse),
            })
        }

        if let Some(unstranded) = counts.unstranded {
            let mismatches = self.mismatches(&sequence, unstranded);
            result.push(IntervalSummary {
                interval: content.interval.clone(),
                strand: self.strandpred.predict(&content.interval, &mismatches),
                name: name.to_owned(),
                mismatches: mismatches,
            })
        }

        result
    }
}

impl<RefNucPred: RefNucPredictor, StrandingPred: LocusStrandPredictor> LocusSummarizer
    for BaseSummarizer<RefNucPred, StrandingPred>
{
    fn summarize(&self, content: NucCounterContent) -> Vec<LocusSummary> {
        // 1. get assembly data
        let reference = self.reference(&content.interval);

        // 2. predict reference nucleotides
        let counts = &content.counts;
        if counts.forward.is_none() && counts.reverse.is_none() && counts.unstranded.is_none() {
            // TODO: Log warning
            return vec![];
        }
        let sequence = self.sequence(&reference, counts);
        let coordinates = &content
            .interval
            .range()
            .map(|x| Locus::new(content.interval.contig().to_owned(), x))
            .collect_vec();

        let mut result = Vec::with_capacity(counts.capacity());
        if let Some(forward) = counts.forward {
            result.extend(izip!(coordinates.clone(), sequence.clone(), forward).map(
                |(locus, refnuc, count)| LocusSummary {
                    locus,
                    strand: Strand::Forward,
                    refnuc,
                    sequenced: *count,
                },
            ))
        }
        if let Some(reverse) = counts.reverse {
            result.extend(izip!(coordinates.clone(), sequence.clone(), reverse).map(
                |(locus, refnuc, count)| LocusSummary {
                    locus,
                    strand: Strand::Reverse,
                    refnuc,
                    sequenced: *count,
                },
            ))
        }

        if let Some(unstranded) = counts.unstranded {
            result.extend(
                izip!(coordinates.clone(), sequence.clone(), unstranded).map(
                    |(locus, refnuc, count)| {
                        let strand = self.strandpred.predict(&locus, &refnuc, count);
                        LocusSummary {
                            locus,
                            strand,
                            refnuc,
                            sequenced: *count,
                        }
                    },
                ),
            )
        }
        result
    }
}
