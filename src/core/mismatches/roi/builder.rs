use std::iter::zip;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::ROIRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::roi::{ROIData, ROIDataVec, ROIMismatchesVec, ROINucCounts};
use crate::core::mismatches::{Batch, Builder};
use crate::core::refpred::{PredNucleotide, RefEngine, RefEngineResult};
use crate::core::rpileup::ncounter::NucCounterResult;
use crate::core::strandutil::Stranded;
use crate::core::workload::ROI;

#[derive(Clone)]
pub struct ROIMismatchesBuilder<RR: ROIRetainer, MP: MismatchesPreFilter<ROIData>> {
    buffer: Vec<NucCounts>,
    refpred: Box<dyn RefEngine>,
    retainer: Option<RR>,
    prefilter: Option<MP>,
}

impl<'a, RR, MP> ROIMismatchesBuilder<RR, MP>
where
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIData>,
{
    pub fn new(maxsize: usize, refpred: Box<dyn RefEngine>, retainer: Option<RR>, prefilter: Option<MP>) -> Self {
        Self { buffer: Vec::with_capacity(maxsize), refpred, retainer, prefilter }
    }

    #[allow(clippy::too_many_arguments)]
    fn process(
        &self,
        cntstart: Position,
        cnts: &'a [NucCounts],
        refpred: &RefEngineResult<'_>,
        roi: &'a ROI,
        coverage: u32,
        retain: &mut ROIDataVec,
        other: &mut ROIDataVec,
    ) {
        // Get mismatches
        let (prednuc, mismatches, heterozygous) = self.summarize(roi, cntstart, refpred.predicted, cnts);
        let record = ROIData { roi: roi.into(), coverage, homozygous: prednuc, heterozygous, mismatches };
        if self.retainer.as_ref().map_or(false, |x| x.retained(roi.contig(), &roi.range(), roi.strand(), roi.name())) {
            // Must be retained
            retain.push(record);
        } else if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&record)) {
            // Must be other
            other.push(record);
        }
    }

    fn summarize(
        &self,
        roi: &'a ROI,
        cntstart: Position,
        prednuc: &'a [PredNucleotide],
        cnts: &'a [NucCounts],
    ) -> (NucCounts, ROINucCounts, u64) {
        debug_assert!(roi.range().start >= cntstart && roi.range().end <= (cntstart + cnts.len() as u64));
        let mut mismatches = ROINucCounts::zeros();
        let mut nuccnts = NucCounts::zeros();
        let mut heterozygous = 0;

        for sub in roi.subintervals() {
            let idx = (sub.start - cntstart) as usize..(sub.end - cntstart) as usize;
            for (nuc, seq) in zip(&prednuc[idx.clone()], &cnts[idx]) {
                match nuc {
                    PredNucleotide::Homozygous(nuc) => match nuc {
                        Nucleotide::A => {
                            nuccnts.A += 1;
                            mismatches.A += seq.into();
                        }
                        Nucleotide::C => {
                            nuccnts.C += 1;
                            mismatches.C += seq.into();
                        }
                        Nucleotide::G => {
                            nuccnts.G += 1;
                            mismatches.G += seq.into();
                        }
                        Nucleotide::T => {
                            nuccnts.T += 1;
                            mismatches.T += seq.into();
                        }
                        // Skip unknown nucleotides
                        Nucleotide::Unknown => {}
                    },
                    PredNucleotide::Heterozygous(_) => {
                        heterozygous += 1;
                        // // Skip sites with unknown alleles
                        // let (n1, n2) = ((*n1).try_into(), (*n2).try_into());
                        // if n1.is_err() || n2.is_err() {
                        //     continue;
                        // }
                        // let (n1, n2): (ReqNucleotide, ReqNucleotide) = (n1.unwrap(), n2.unwrap());
                        // // Fill in heterozygous matches / mismatches
                        // let matches = std::cmp::min(seq[n1], seq[n2]) as f32;
                        // mismatches[n1][n1] += matches;
                        // mismatches[n2][n2] += matches;
                        // // Unbalanced alleles counts
                        // if seq[n1] > seq[n2] {
                        //     mismatches[n2][n1] += (seq[n1] - seq[n2]) as f32;
                        // } else {
                        //     mismatches[n1][n2] += (seq[n2] - seq[n1]) as f32;
                        // }
                        // // Half-weight other mismatches
                        // for n in [n1, n2] {
                        //     for subs in [ReqNucleotide::A, ReqNucleotide::C, ReqNucleotide::G, ReqNucleotide::T] {
                        //         if subs == n1 || subs == n2 || seq[subs] == 0 {
                        //             continue;
                        //         }
                        //         mismatches[n][subs] += (seq[subs] as f32) * 0.5;
                        //     }
                        // }
                    }
                }
            }
        }

        (nuccnts, mismatches, heterozygous)
    }

    #[inline]
    fn size_hint(&self, nc: &NucCounterResult<'a, &'a ROI>) -> Stranded<usize> {
        let mut size: Stranded<usize> = Stranded::default();
        for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
            size[strand] = nc.cnts.iter().map(|item| item.cnts[strand].is_some() as usize).sum()
        }
        size
    }
}

impl<'a, RR, MP> Builder<'a> for ROIMismatchesBuilder<RR, MP>
where
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIData>,
{
    type Out = ROIMismatchesVec;
    type SourceCounts = NucCounterResult<'a, &'a ROI>;

    fn build(&mut self, nc: Self::SourceCounts) -> Batch<Self::Out> {
        let contig = nc.contig.to_owned();

        // Pre-allocate results
        let _hint = self.size_hint(&nc);
        let mut items = Stranded::with_fn(|strnd| {
            ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::new())
            // ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::with_capacity(hint[strnd]))
        });
        let mut retained = Stranded::with_fn(|strnd| {
            ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::new())
            // ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::with_capacity(hint[strnd] / 10))
        });

        for item in nc.cnts.into_iter() {
            // debug_assert!(item.coverage.forward + item.coverage.reverse + item.coverage.unknown > 0);

            // Predict the reference
            let counts = item.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(&contig, item.range.clone(), counts);
            let refpred = self.refpred.results();

            // Process the counts
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                if let Some(cnts) = item.cnts[strand] {
                    self.process(
                        item.range.start,
                        cnts,
                        &refpred,
                        item.data,
                        item.coverage[strand],
                        &mut retained[strand].data,
                        &mut items[strand].data,
                    );
                }
            }
        }
        Batch { contig, mapped: nc.mapped, retained, items }
    }
}
