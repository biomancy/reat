use std::path::PathBuf;

use bio_types::genome::AbstractInterval;

use itertools::zip;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::core::counting::nuccounter::NucCounter;

use crate::core::dna::Nucleotide;

use crate::core::refnuc::RefNucPredictor;
use crate::core::runner::{RunResult, Runner};
use crate::core::workload::ROIWorkload;

// Dummy FastaReader that is safe to Send
pub struct FastaReader {
    pub faidx: faidx::Reader,
    pub path: PathBuf,
}
unsafe impl Send for FastaReader {}

pub struct BaseRunner<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor> {
    htsreaders: Vec<bam::IndexedReader>,
    htsfiles: Vec<PathBuf>,
    refreader: FastaReader,
    mapped: u32,
    counter: Counter,
    refnucpred: RefNucPred,
    cachenuc: Vec<Nucleotide>,
}

impl<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor> BaseRunner<Counter, RefNucPred> {
    pub fn new(htsfiles: Vec<PathBuf>, reference: PathBuf, counter: Counter, refnucpred: RefNucPred) -> Self {
        let htsreaders: Vec<bam::IndexedReader> = htsfiles
            .iter()
            .map(|hts| {
                bam::IndexedReader::from_path(&hts).unwrap_or_else(|_| panic!("Failed to open file {}", hts.display()))
            })
            .collect();
        let refreader = faidx::Reader::from_path(reference.as_path())
            .unwrap_or_else(|_| panic!("Failed to open file {}", reference.display()));
        Self {
            htsreaders,
            htsfiles,
            refreader: FastaReader { faidx: refreader, path: reference },
            mapped: 0,
            counter,
            refnucpred,
            cachenuc: Vec::new(),
        }
    }

    fn predseq(&mut self) {
        self.cachenuc.clear();

        let interval = self.counter.interval();
        let (start, end) = (interval.range().start as usize, interval.range().end as usize);
        let reference = self
            .refreader
            .faidx
            .fetch_seq(interval.contig(), start, end - 1)
            .unwrap_or_else(|_| panic!("Failed to fetch sequence for region {:?}", interval));

        debug_assert_eq!(reference.len(), end - start);
        let counted = self.counter.counted();
        debug_assert_eq!(reference.len(), counted.len());
        let refpred = &self.refnucpred;
        let iter = zip(reference, counted)
            .map(|(a, b)| (Nucleotide::from(*a), b))
            .map(|(refnuc, seqnuc)| refpred.predict(&refnuc, seqnuc));
        self.cachenuc.extend(iter);
    }

    fn reset(&mut self, workload: ROIWorkload) {
        self.counter.reset(workload);
        self.cachenuc.clear();
    }

    fn count(&mut self) {
        for htsreader in &mut self.htsreaders {
            let interval = self.counter.interval();
            if !htsreader.header().target_names().contains(&interval.contig().as_bytes()) {
                continue;
            }

            htsreader
                .fetch((interval.contig(), interval.range().start, interval.range().end))
                .unwrap_or_else(|_| panic!("Failed to fetch reads for {:?} (HTS file corrupted?)", interval));

            let mut record = Record::new();
            while let Some(r) = htsreader.read(&mut record) {
                r.expect("Failed to parse record information (HTS file corrupted?)");
                self.counter.count(&mut record);
            }
        }
    }
}

impl<Counter: NucCounter<Record> + Clone, RefNucPred: RefNucPredictor + Clone> Clone
    for BaseRunner<Counter, RefNucPred>
{
    fn clone(&self) -> Self {
        Self::new(self.htsfiles.clone(), self.refreader.path.clone(), self.counter.clone(), self.refnucpred.clone())
    }
}

impl<Counter: NucCounter<Record>, RefNucPred: RefNucPredictor> Runner for BaseRunner<Counter, RefNucPred> {
    fn run(&mut self, workload: ROIWorkload) -> Vec<RunResult> {
        self.reset(workload);
        self.count();
        self.mapped += self.counter.mapped();
        self.predseq();

        let sequence = &self.cachenuc;
        let counted = self.counter.results();
        let interstart = self.counter.interval().range().start;

        let mut results = Vec::with_capacity(counted.len() * 2);
        for c in counted {
            let (start, end) = (c.roi.range().start - interstart, c.roi.range().end - interstart);
            let (start, end) = (start as usize, end as usize);
            let runres = RunResult {
                interval: c.roi,
                name: c.name,
                strand: c.strand,
                cnts: c.cnts,
                reference: &sequence[start..end],
            };
            debug_assert_eq!(runres.reference.len(), runres.cnts.nuc.len());
            results.push(runres);
        }
        results
    }

    #[inline]
    fn mapped(&self) -> u32 {
        self.mapped
    }
}
