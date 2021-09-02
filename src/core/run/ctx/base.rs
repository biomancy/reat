use std::path::{Path, PathBuf};

use bio_types::genome::{AbstractInterval, Interval};
use itertools::izip;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::core::counting::{CountsBufferContent, NucCounter, NucCounterContent, NucCounts};
use crate::core::dna::Nucleotide;
use crate::core::filtering::summary::{IntervalSummaryFilter, LocusSummaryFilter};
use crate::core::refnuc::RefNucPredictor;
use crate::core::run::ctx::{LociRunCtx, ROIRunCtx, RunCtx};
use crate::core::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};

// Dummy FastaReader that is safe to Send
pub struct FastaReader(faidx::Reader, PathBuf);
unsafe impl Send for FastaReader {}

pub struct BaseRunCtx<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor, StrandPred, Filter> {
    htsreaders: Vec<bam::IndexedReader>,
    htsfiles: Vec<PathBuf>,
    refreader: FastaReader,
    reads_counted: u32,
    counter: Counter,
    refnucpred: RefNucPred,
    strandpred: StrandPred,
    filter: Filter,
}

impl<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor, StrandPred, Filter>
    BaseRunCtx<Counter, RefNucPred, StrandPred, Filter>
{
    pub fn new(
        htsfiles: &[PathBuf],
        reference: &Path,
        counter: Counter,
        refnucpred: RefNucPred,
        strandpred: StrandPred,
        filter: Filter,
    ) -> Self {
        let htsreaders: Vec<bam::IndexedReader> = htsfiles
            .iter()
            .map(|hts| {
                bam::IndexedReader::from_path(&hts).unwrap_or_else(|_| panic!("Failed to open file {}", hts.display()))
            })
            .collect();
        let refreader = faidx::Reader::from_path(reference)
            .unwrap_or_else(|_| panic!("Failed to open file {}", reference.display()));
        BaseRunCtx {
            htsreaders,
            htsfiles: htsfiles.to_vec(),
            refreader: FastaReader(refreader, reference.into()),
            reads_counted: 0,
            counter,
            refnucpred,
            strandpred,
            filter,
        }
    }

    fn implpredseq(
        refnucpred: &impl RefNucPredictor,
        reference: &[Nucleotide],
        counts: &CountsBufferContent,
    ) -> Vec<Nucleotide> {
        let sumup = |a: &[NucCounts], b: &[NucCounts]| izip!(a, b).map(|(a, b)| *a + *b).collect();

        let storage: Vec<NucCounts>;
        let sequenced: &[NucCounts] = match (counts.forward, counts.reverse, counts.unstranded) {
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

    pub fn assembly(&self, interval: &Interval) -> Vec<Nucleotide> {
        let (start, end) = (interval.range().start as usize, interval.range().end as usize);
        let reference = self
            .refreader
            .0
            .fetch_seq(interval.contig(), start, end - 1)
            .unwrap_or_else(|_| panic!("Failed to fetch sequence for region {:?}", interval));

        debug_assert_eq!(reference.len(), end - start);

        reference.iter().map(|e| Nucleotide::from(*e)).collect()
    }
}

impl<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor, StrandPred, Filter> RunCtx
    for BaseRunCtx<Counter, RefNucPred, StrandPred, Filter>
{
    fn count(&mut self, interval: &Interval) {
        self.counter.reset(interval.clone());

        for htsreader in &mut self.htsreaders {
            if !htsreader.header().target_names().contains(&interval.contig().as_bytes()) {
                continue;
            }

            htsreader
                .fetch((interval.contig(), interval.range().start, interval.range().end))
                .unwrap_or_else(|_| panic!("Failed to fetch reads for {:?} (HTS file corrupted?)", interval));

            let mut record = Record::new();
            while let Some(r) = htsreader.read(&mut record) {
                r.expect("Failed to parse record information (HTS file corrupted?)");
                self.counter.process(&mut record);
            }
        }
        self.reads_counted += self.counter.reads_counted()
    }

    fn finalize(&self) -> Option<(NucCounterContent, Vec<Nucleotide>)> {
        if self.counter.empty() {
            None
        } else {
            let nuccontent = self.counter.content();
            let assembly = self.assembly(&nuccontent.interval);
            let sequence = self.predseq(&assembly, &nuccontent.counts);
            Some((nuccontent, sequence))
        }
    }

    fn htsfiles(&self) -> &[PathBuf] {
        &self.htsfiles
    }

    fn reads_counted(&self) -> u32 {
        self.reads_counted
    }
}

impl<Counter, RefNucPred, StrandPred, Filter> Clone for BaseRunCtx<Counter, RefNucPred, StrandPred, Filter>
where
    Counter: NucCounter<Record> + Clone,
    RefNucPred: RefNucPredictor + Clone,
    StrandPred: Clone,
    Filter: Clone,
{
    fn clone(&self) -> Self {
        Self::new(
            &self.htsfiles,
            &self.refreader.1,
            self.counter.clone(),
            self.refnucpred.clone(),
            self.strandpred.clone(),
            self.filter.clone(),
        )
    }
}

impl<Counter, RefNucPred, StrandPred, Filter> ROIRunCtx for BaseRunCtx<Counter, RefNucPred, StrandPred, Filter>
where
    Counter: NucCounter<Record> + Clone,
    RefNucPred: RefNucPredictor + Clone,
    StrandPred: IntervalStrandPredictor,
    Filter: IntervalSummaryFilter,
{
    type StrandPred = StrandPred;
    type Filter = Filter;

    #[inline]
    fn strandpred(&self) -> &Self::StrandPred {
        &self.strandpred
    }

    #[inline]
    fn filter(&self) -> &Self::Filter {
        &self.filter
    }
}

impl<Counter, RefNucPred, StrandPred, Filter> LociRunCtx for BaseRunCtx<Counter, RefNucPred, StrandPred, Filter>
where
    Counter: NucCounter<Record> + Clone,
    RefNucPred: RefNucPredictor + Clone,
    StrandPred: LocusStrandPredictor,
    Filter: LocusSummaryFilter,
{
    type StrandPred = StrandPred;
    type Filter = Filter;

    #[inline]
    fn strandpred(&self) -> &Self::StrandPred {
        &self.strandpred
    }

    #[inline]
    fn filter(&self) -> &Self::Filter {
        &self.filter
    }
}

#[cfg(test)]
mod tests {
    use mockall::predicate::eq;

    use crate::core::counting::MockNucCounter;
    use crate::core::filtering::summary::MockIntervalSummaryFilter;
    use crate::core::refnuc::MockRefNucPredictor;
    use crate::core::stranding::predict::MockIntervalStrandPredictor;

    use super::*;

    type DummyContext =
        BaseRunCtx<MockNucCounter<Record>, MockRefNucPredictor, MockIntervalStrandPredictor, MockIntervalSummaryFilter>;

    #[test]
    fn predsequence() {
        let mut refnucpred = MockRefNucPredictor::new();
        let reference = vec![Nucleotide::A, Nucleotide::C, Nucleotide::G];
        let counts = NucCounts { A: 1, C: 2, G: 3, T: 4 };
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
            let pred = DummyContext::implpredseq(&refnucpred, &reference, &content);
            assert_eq!(pred, reference);
        }

        let content = CountsBufferContent { forward: None, reverse: None, unstranded: None };
        let pred = DummyContext::implpredseq(&refnucpred, &reference, &content);
        assert_eq!(pred, vec![]);
    }
}
