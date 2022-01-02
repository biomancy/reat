use std::marker::PhantomData;
use std::path::PathBuf;

use bio_types::genome::AbstractInterval;
use bio_types::genome::Interval;
use derive_more::Constructor;
use itertools::Itertools;
use rust_htslib::bam;
use rust_htslib::bam::pileup::Pileup;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

use crate::core::dna::Nucleotide;
use crate::core::hooks::ProcessingHook;
use crate::core::mismatches::IntermediateMismatches;
use crate::core::refpred::RefEngine;
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::{NucCounter, NucCountingResult};
use crate::core::rpileup::{ReadsCollider, ReadsPileupEngine};
use crate::core::stranding::predict::StrandingEngine;

pub struct REATRunner<'a, Summary, NCounter, RefPred, Strander>
where
    Summary: IntermediateMismatches,
    NCounter: NucCounter<'a, Record, Summary>,
    RefPred: RefEngine,
    Strander: StrandingEngine<Summary>,
    <NCounter as ReadsCollider<'a, Record>>::ColliderResult: NucCountingResult<'a, Summary>,
{
    pileuper: HTSPileupEngine<'a, NCounter>,
    refpred: RefPred,
    strander: Strander,
    hooks: Vec<Box<dyn ProcessingHook<Summary>>>,
    phantom: PhantomData<&'a u8>,
}

impl<'a, Summary, NCounter, RefPred, Strander> REATRunner<'a, Summary, NCounter, RefPred, Strander>
where
    Summary: IntermediateMismatches,
    NCounter: NucCounter<'a, Record, Summary>,
    RefPred: RefEngine,
    Strander: StrandingEngine<Summary>,
    <NCounter as ReadsCollider<'a, Record>>::ColliderResult: NucCountingResult<'a, Summary>,
{
    pub fn new(
        refpred: RefPred,
        strander: Strander,
        hooks: Vec<Box<dyn ProcessingHook<Summary>>>,
        pileuper: HTSPileupEngine<'a, NCounter>,
    ) -> Self {
        Self { pileuper, refpred, strander, hooks, phantom: Default::default() }
    }

    pub fn run(&'a mut self, interval: Interval, workload: NCounter::Workload) -> Vec<Summary::Finalized> {
        self.pileuper.run(interval, workload);
        let result = self.pileuper.result().unwrap();

        self.refpred.reset();
        for r in &result {
            self.refpred.run(r.interval(), r.ncounts());
        }
        let reference = self.refpred.results();

        let result = result.into_iter().zip(reference).map(|(x, seq)| x.mismatches(seq)).flatten().collect_vec();

        // Run stranding
        let mut result = self.strander.strand(result);

        // Pre-final hooks
        for x in &mut self.hooks {
            result = x.hook(result);
        }

        // Finalize and return
        result.into_iter().map(|x| x.finalize()).collect()
    }
}
