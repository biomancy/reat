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
use crate::core::hooks::HooksEngine;
use crate::core::mismatches::IntermediateMismatches;
use crate::core::refpred::RefEngine;
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::{CountingResults, NucCounter, ToMismatches};
use crate::core::rpileup::{ReadsCollider, ReadsPileupEngine};
use crate::core::stranding::predict::StrandingEngine;

pub struct REATRunner<'a, Mismatches, RefPred, Strander, Hook, NCounter>
where
    Mismatches: IntermediateMismatches,
    RefPred: RefEngine,
    Strander: StrandingEngine<Mismatches>,
    Hook: HooksEngine<Mismatches>,
    NCounter: NucCounter<'a, Record>,
    NCounter::NucCounts: ToMismatches<'a, Mismatches>,
{
    pileuper: HTSPileupEngine<'a, NCounter>,
    refpred: RefPred,
    strander: Strander,
    hooks: Hook,
    phantom: PhantomData<&'a Mismatches>,
}

impl<'a, Mismatches, RefPred, Strander, Hook, NCounter> REATRunner<'a, Mismatches, RefPred, Strander, Hook, NCounter>
where
    Mismatches: IntermediateMismatches,
    RefPred: RefEngine,
    Strander: StrandingEngine<Mismatches>,
    Hook: HooksEngine<Mismatches>,
    NCounter: NucCounter<'a, Record>,
    NCounter::NucCounts: ToMismatches<'a, Mismatches>,
{
    pub fn new(refpred: RefPred, strander: Strander, hooks: Hook, pileuper: HTSPileupEngine<'a, NCounter>) -> Self {
        Self { pileuper, refpred, strander, hooks, phantom: Default::default() }
    }

    pub fn run(&'a mut self, workload: NCounter::Workload) -> Vec<Mismatches::Finalized> {
        self.pileuper.run(workload);
        let result = self.pileuper.result().unwrap();

        self.refpred.reset();
        for r in result.items() {
            let interval = Interval::new(r.contig().into(), r.range().clone());
            self.refpred.run(&interval, r.ncounts());
        }
        let reference = self.refpred.results();

        let result =
            result.dissolve().1.into_iter().zip(reference).map(|(x, seq)| x.mismatches(seq)).flatten().collect_vec();

        // Run stranding
        let mut result = self.strander.strand(result);

        // Final hooks
        result = self.hooks.on_finish(result);

        // Finalize and return
        result.into_iter().map(|x| x.finalize()).collect()
    }
}
