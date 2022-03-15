

use rust_htslib::bam::Record;

use crate::core::hooks::stats::EditingStat;
use crate::core::hooks::HooksEngine;
use crate::core::mismatches::{Batch, Builder, MismatchesVec};
use crate::core::rpileup::hts::HTSPileupEngine;

use crate::core::rpileup::{ReadsCollider, ReadsCollidingEngine};
use crate::core::stranding::predict::StrandingEngine;

pub trait Runner<'runner, T: MismatchesVec> {
    type Workload;

    fn run(&'runner mut self, workload: Self::Workload) -> Option<Batch<T>>;
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

impl<'runner, NCounter, MBuilder, Strander, Hook> Runner<'runner, MBuilder::Out>
    for REATRunner<NCounter, MBuilder, Strander, Hook>
where
    for<'a> NCounter: ReadsCollider<'a, Record>,
    MBuilder: Builder<'runner, SourceCounts = <NCounter as ReadsCollider<'runner, Record>>::ColliderResult>,
    Strander: StrandingEngine<MBuilder::Out>,
    Hook: HooksEngine<MBuilder::Out> + Clone,
{
    type Workload = <NCounter as ReadsCollider<'runner, Record>>::Workload;

    fn run(&'runner mut self, workload: Self::Workload) -> Option<Batch<MBuilder::Out>> {
        self.pileuper.run(workload);

        let ncounts = match self.pileuper.result() {
            Some(x) => x,
            None => return None,
        };

        let mut batch = self.mmbuilder.build(ncounts);

        // Run stranding
        batch.retained = self.strander.strand(&batch.contig, batch.retained);
        batch.items = self.strander.strand(&batch.contig, batch.items);

        // Final hooks
        self.hook.on_finish(&mut batch);
        Some(batch)
    }

    fn stats(self) -> Vec<Box<dyn EditingStat<MBuilder::Out>>> {
        self.hook.stats()
    }
}
