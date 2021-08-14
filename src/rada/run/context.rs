use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bio_types::genome::{AbstractInterval, Interval};
use itertools::izip;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::rada::counting::{CountsBufferContent, LocusCounts, NucCounter};
use crate::rada::dna::Nucleotide;
use crate::rada::read::AlignedRead;
use crate::rada::refnuc::RefNucPredictor;

pub struct FastaReader(faidx::Reader);
unsafe impl Send for FastaReader {}

pub struct ThreadContext<R: AlignedRead, Counter: NucCounter<R>, RefNucPred: RefNucPredictor, StrandPred, SummaryFilter>
{
    pub(super) htsreaders: Vec<bam::IndexedReader>,
    pub(super) reference: FastaReader,
    pub(super) counter: Counter,
    pub(super) refnucpred: RefNucPred,
    pub(super) strandpred: StrandPred,
    pub(super) filter: SummaryFilter,
    pub(super) phantom: PhantomData<fn() -> R>,
}

impl<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor, StrandPred, SummaryFilter>
    ThreadContext<Record, Counter, RefNucPred, StrandPred, SummaryFilter>
{
    pub fn new(
        htsfiles: &[PathBuf],
        reference: &Path,
        counter: Counter,
        refnucpred: RefNucPred,
        strandpred: StrandPred,
        filter: SummaryFilter,
    ) -> Self {
        let htsreaders: Vec<bam::IndexedReader> = htsfiles
            .iter()
            .map(|hts| {
                bam::IndexedReader::from_path(&hts).unwrap_or_else(|_| panic!("Failed to open file {}", hts.display()))
            })
            .collect();
        let reference = faidx::Reader::from_path(reference)
            .unwrap_or_else(|_| panic!("Failed to open file {}", reference.display()));
        ThreadContext {
            htsreaders,
            reference: FastaReader(reference),
            counter,
            refnucpred,
            strandpred,
            filter,
            phantom: Default::default(),
        }
    }

    pub fn nuccount(&mut self, interval: &Interval) {
        self.counter.reset(interval.clone());

        for htsreader in &mut self.htsreaders {
            let err = htsreader.fetch((interval.contig(), interval.range().start, interval.range().end));
            err.unwrap_or_else(|_| panic!("Failed to fetch reads for {:?} (HTS file corrupted?)", interval));

            let mut record = Record::new();
            while let Some(r) = htsreader.read(&mut record) {
                r.expect("Failed to parse record information (HTS file corrupted?)");
                self.counter.process(&mut record);
            }
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
        let (start, end) = (interval.range().start as usize, interval.range().end as usize);
        let reference = self
            .reference
            .0
            .fetch_seq(interval.contig(), start, end - 1)
            .unwrap_or_else(|_| panic!("Failed to fetch sequence for region {:?}", interval));

        debug_assert_eq!(reference.len(), end - start);

        reference.iter().map(|e| Nucleotide::from(*e)).collect()
    }
}

#[cfg(test)]
mod tests {
    use mockall::predicate::eq;

    use crate::rada::counting::MockNucCounter;
    use crate::rada::filtering::summary::MockIntervalSummaryFilter;
    use crate::rada::refnuc::MockRefNucPredictor;
    use crate::rada::stranding::predict::MockIntervalStrandPredictor;

    use super::*;

    type Context = ThreadContext<
        Record,
        MockNucCounter<Record>,
        MockRefNucPredictor,
        MockIntervalStrandPredictor,
        MockIntervalSummaryFilter,
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
