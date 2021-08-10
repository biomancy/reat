use std::collections::HashMap;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{Same, Strand};
use itertools::izip;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::rada::counting::{CountsBufferContent, LocusCounts, NucCounter, NucCounterContent};
use crate::rada::dna::Nucleotide;
use crate::rada::filtering::summary::{IntervalSummaryFilter, LocusSummaryFilter};
use crate::rada::read::AlignedRead;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor, StrandPredictor};
use crate::rada::summary::IntervalSummary;
use crate::rada::workload::Workload;

pub struct ThreadContext<R: AlignedRead, SummaryFilter, Counter: NucCounter<R>, RefNucPred: RefNucPredictor, StrandPred>
{
    pub(super) htsreaders: HashMap<PathBuf, bam::IndexedReader>,
    pub(super) reference: faidx::Reader,
    pub(super) filter: SummaryFilter,
    pub(super) counter: Counter,
    pub(super) refnucpred: RefNucPred,
    pub(super) strandpred: StrandPred,
    pub(super) phantom: PhantomData<R>,
}

impl<SummaryFilter, Counter: NucCounter<Record>, RefNucPred: RefNucPredictor, StrandPred>
    ThreadContext<Record, SummaryFilter, Counter, RefNucPred, StrandPred>
{
    pub fn new(
        htsfiles: &[&Path],
        reference: &Path,
        filter: SummaryFilter,
        counter: Counter,
        refnucpred: RefNucPred,
        strandpred: StrandPred,
    ) -> Self {
        let mut htsreaders: HashMap<PathBuf, bam::IndexedReader> = HashMap::new();
        for &hts in htsfiles {
            let reader = bam::IndexedReader::from_path(&hts).expect(&format!("Failed to open file {}", hts.display()));
            htsreaders.insert(PathBuf::from(hts), reader);
        }
        let reference =
            faidx::Reader::from_path(reference).expect(&format!("Failed to open file {}", reference.display()));
        ThreadContext { htsreaders, reference, filter, counter, refnucpred, strandpred, phantom: Default::default() }
    }

    pub fn nuccount(&mut self, bam: &Path, interval: &Interval) {
        self.counter.reset(interval.clone());

        let htsreader = self.htsreaders.get_mut(bam).expect("Bug: thread failed to initialize all BAM readers");

        let err = htsreader.fetch((interval.contig(), interval.range().start, interval.range().end));
        err.expect(&format!("Failed to fetch reads for {:?} (HTS file corrupted?)", interval));

        let mut record = Record::new();
        while let Some(r) = htsreader.read(&mut record) {
            r.expect("Failed to parse record information (HTS file corrupted?)");
            self.counter.process(&mut record);
        }
    }

    fn implpredseq(
        refnucpred: &impl RefNucPredictor,
        reference: &[Nucleotide],
        counts: &CountsBufferContent,
    ) -> Vec<Nucleotide> {
        let sumup = |a: &[LocusCounts], b: &[LocusCounts]| izip!(a, b).map(|(a, b)| *a + *b).collect();

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
                storage = izip!(f, r, u).map(|(a, b, c)| *a + *b + *c).collect();
                &storage
            }
            (None, None, None) => return vec![],
        };

        izip!(reference, sequenced).map(|(refnuc, seqnuc)| refnucpred.predict(refnuc, seqnuc)).collect()
    }

    pub fn predseq(&self, reference: &[Nucleotide], counts: &CountsBufferContent) -> Vec<Nucleotide> {
        Self::implpredseq(&self.refnucpred, reference, counts)
    }

    pub fn reference(&self, interval: &Interval) -> Vec<Nucleotide> {
        let (start, end) = (interval.range().start as usize, interval.range().end as usize - 1);
        let reference = self
            .reference
            .fetch_seq(interval.contig(), start, end)
            .expect(&format!("Failed to fetch sequence for region {:?}", interval));

        debug_assert_eq!(reference.len(), end - start);

        reference.iter().map(|e| Nucleotide::from(*e)).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::rada::counting::MockNucCounter;
    use crate::rada::filtering::summary::MockIntervalSummaryFilter;
    use crate::rada::refnuc::MockRefNucPredictor;
    use crate::rada::stranding::predict::MockIntervalStrandPredictor;
    use bio_types::strand::Same;
    use mockall::predicate::eq;
    use rust_htslib::errors::Error::ThreadPool;

    use super::*;

    type Context = ThreadContext<
        Record,
        MockIntervalSummaryFilter,
        MockNucCounter<Record>,
        MockRefNucPredictor,
        MockIntervalStrandPredictor,
    >;

    #[test]
    fn predsequence() {
        let mut refnucpred = MockRefNucPredictor::new();
        let reference = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];
        let counts = LocusCounts { A: 1, C: 2, G: 3, T: 4 };
        let sequenced = vec![counts].repeat(3);
        let content = Some(sequenced.as_slice());

        for (forward, reverse, unstranded, mult) in [
            (content, content, content, 3u32),
            (content, content, None, 2u32),
            (content, None, content, 2u32),
            (None, content, content, 2u32),
            (content, None, None, 1u32),
            (None, content, None, 1u32),
            (None, None, content, 1u32),
        ] {
            for nuc in &reference {
                refnucpred.expect_predict().once().with(eq(*nuc), eq(counts * mult)).return_const(*nuc);
            }

            let content = CountsBufferContent { forward, reverse, unstranded };
            let pred = Context::implpredseq(&refnucpred, &reference, &content);
            assert_eq!(pred, reference);
        }

        let content = CountsBufferContent { forward: None, reverse: None, unstranded: None };
        let pred = Context::implpredseq(&refnucpred, &reference, &content);
        assert_eq!(pred, vec![]);
    }
}
