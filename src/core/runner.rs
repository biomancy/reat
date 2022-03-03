use itertools::chain;
use rust_htslib::bam::Record;

use crate::core::hooks::stats::EditingStat;
use crate::core::hooks::HooksEngine;
use crate::core::mismatches::{BatchedMismatches, BatchedMismatchesBuilder};
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::rpileup::{ReadsCollider, ReadsCollidingEngine};
use crate::core::stranding::predict::StrandingEngine;

pub trait Runner<'runner, T: BatchedMismatches> {
    type Workload;
    type Result;

    fn run(&'runner mut self, workload: Self::Workload) -> Self::Result;
    fn stats(self) -> Vec<Box<dyn EditingStat<T>>>;
}

#[derive(Clone)]
pub struct REATRunner<NCounter, MismatchesBuilder, Strander, Hook>
where
    for<'a> NCounter: ReadsCollider<'a, Record>,
{
    pileuper: HTSPileupEngine<NCounter>,
    mmbuilder: MismatchesBuilder,
    strander: Strander,
    hook: Hook,
}

impl<NCounter, MismatchesBuilder, Strander, Hook> REATRunner<NCounter, MismatchesBuilder, Strander, Hook>
where
    for<'a> NCounter: ReadsCollider<'a, Record>,
{
    pub fn new(
        mmbuilder: MismatchesBuilder,
        strander: Strander,
        pileuper: HTSPileupEngine<NCounter>,
        hook: Hook,
    ) -> Self {
        Self { pileuper, mmbuilder, strander, hook }
    }
}

impl<'runner, NCounter, MismatchesBuilder, Strander, Hook> Runner<'runner, MismatchesBuilder::Mismatches>
    for REATRunner<NCounter, MismatchesBuilder, Strander, Hook>
where
    for<'a> NCounter: ReadsCollider<'a, Record>,
    <NCounter as ReadsCollider<'runner, Record>>::ColliderResult: AggregatedNucCounts<'runner>,
    MismatchesBuilder: BatchedMismatchesBuilder<'runner, <NCounter as ReadsCollider<'runner, Record>>::ColliderResult>,
    Strander: StrandingEngine<MismatchesBuilder::Mismatches>,
    Hook: HooksEngine<MismatchesBuilder::Mismatches> + Clone,
{
    type Workload = <NCounter as ReadsCollider<'runner, Record>>::Workload;
    type Result = Vec<<MismatchesBuilder::Mismatches as BatchedMismatches>::Flattened>;

    fn run(&'runner mut self, workload: Self::Workload) -> Self::Result {
        self.pileuper.run(workload);

        let ncounts = self.pileuper.result().unwrap();

        let mut mismatches = self.mmbuilder.build(ncounts);

        // Run stranding
        mismatches.retained = self.strander.strand(mismatches.retained);
        mismatches.other = self.strander.strand(mismatches.other);

        // Final hooks
        self.hook.on_finish(&mut mismatches);

        // Finalize and return
        let mismatches = chain(mismatches.retained.into_iter(), mismatches.other.into_iter())
            .map(|x| x.flatten())
            .flatten()
            .collect();
        mismatches
    }

    fn stats(self) -> Vec<Box<dyn EditingStat<MismatchesBuilder::Mismatches>>> {
        self.hook.stats()
    }
}
