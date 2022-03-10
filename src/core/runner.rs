use bio_types::strand::Strand;
use itertools::{chain, Itertools};
use rust_htslib::bam::Record;

use crate::core::hooks::stats::EditingStat;
use crate::core::hooks::HooksEngine;
use crate::core::mismatches::{MismatchesBuilder, MismatchesVec};
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::rpileup::{ReadsCollider, ReadsCollidingEngine};
use crate::core::stranding::predict::StrandingEngine;

pub trait Runner<'runner, T: MismatchesVec> {
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

impl<'runner, NCounter, MBuilder, Strander, Hook> Runner<'runner, MBuilder::Mismatches>
    for REATRunner<NCounter, MBuilder, Strander, Hook>
where
    for<'a> NCounter: ReadsCollider<'a, Record>,
    <NCounter as ReadsCollider<'runner, Record>>::ColliderResult: AggregatedNucCounts<'runner>,
    MBuilder: MismatchesBuilder<'runner, <NCounter as ReadsCollider<'runner, Record>>::ColliderResult>,
    Strander: StrandingEngine<MBuilder::Mismatches>,
    Hook: HooksEngine<MBuilder::Mismatches> + Clone,
{
    type Workload = <NCounter as ReadsCollider<'runner, Record>>::Workload;
    type Result = Vec<<MBuilder::Mismatches as MismatchesVec>::Flat>;

    fn run(&'runner mut self, workload: Self::Workload) -> Self::Result {
        self.pileuper.run(workload);

        let ncounts = match self.pileuper.result() {
            Some(x) => x,
            None => return vec![],
        };

        let mut mismatches = self.mmbuilder.build(ncounts);

        // Run stranding
        // mismatches.retained = self.strander.strand(mismatches.retained);
        // mismatches.items = self.strander.strand(mismatches.items);
        //
        // // Final hooks
        // self.hook.on_finish(&mut mismatches);
        //
        // // Finalize and return
        // chain!(
        //     mismatches.retained.forward.into_iter(),
        //     mismatches.retained.reverse.into_iter(),
        //     mismatches.retained.unknown.into_iter(),
        //     mismatches.items.forward.into_iter(),
        //     mismatches.items.reverse.into_iter(),
        //     mismatches.items.unknown.into_iter(),
        // )
        // .map(|x| x.flatten())
        // .flatten()
        // .collect()
        todo!()
    }

    fn stats(self) -> Vec<Box<dyn EditingStat<MBuilder::Mismatches>>> {
        self.hook.stats()
    }
}
