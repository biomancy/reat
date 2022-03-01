use std::marker::PhantomData;
use std::path::PathBuf;

use bio_types::genome::AbstractInterval;
use bio_types::genome::Interval;
use derive_more::Constructor;
use itertools::{chain, Itertools};
use rust_htslib::bam;
use rust_htslib::bam::pileup::Pileup;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

use crate::core::dna::Nucleotide;
use crate::core::hooks::HooksEngine;
use crate::core::mismatches::{BatchedMismatches, BatchedMismatchesBuilder};
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, NucCounter};
use crate::core::rpileup::ReadsCollidingEngine;
use crate::core::stranding::predict::StrandingEngine;

pub struct REATRunner<'a, NCounter, MismatchesBuilder, Strander, Hook>
where
    NCounter: NucCounter<'a, Record>,
    NCounter::ColliderResult: AggregatedNucCounts<'a>,
    MismatchesBuilder: BatchedMismatchesBuilder<'a, NCounter::ColliderResult>,
    Strander: StrandingEngine<MismatchesBuilder::Mismatches>,
    Hook: HooksEngine<MismatchesBuilder::Mismatches>,
{
    pileuper: HTSPileupEngine<'a, NCounter>,
    mmbuilder: MismatchesBuilder,
    strander: Strander,
    hook: Hook,
}

impl<'a, NCounter, MismatchesBuilder, Strander, Hook> REATRunner<'a, NCounter, MismatchesBuilder, Strander, Hook>
where
    NCounter: NucCounter<'a, Record>,
    NCounter::ColliderResult: AggregatedNucCounts<'a>,
    MismatchesBuilder: BatchedMismatchesBuilder<'a, NCounter::ColliderResult>,
    Strander: StrandingEngine<MismatchesBuilder::Mismatches>,
    Hook: HooksEngine<MismatchesBuilder::Mismatches>,
{
    pub fn new(
        mmbuilder: MismatchesBuilder,
        strander: Strander,
        pileuper: HTSPileupEngine<'a, NCounter>,
        hook: Hook,
    ) -> Self {
        Self { pileuper, mmbuilder, strander, hook }
    }

    pub fn run(
        &'a mut self,
        workload: NCounter::Workload,
    ) -> Vec<<MismatchesBuilder::Mismatches as BatchedMismatches>::Flattened> {
        self.pileuper.run(workload);
        let ncounts = self.pileuper.result().unwrap();

        let mut mismatches = self.mmbuilder.build(ncounts);

        self.hook.on_created(&mut mismatches);

        // Run stranding
        mismatches.retained = self.strander.strand(mismatches.retained);
        mismatches.other = self.strander.strand(mismatches.other);

        self.hook.on_stranded(&mut mismatches);

        // Final hooks
        self.hook.on_finish(&mut mismatches);

        // Finalize and return
        chain(mismatches.retained.into_iter(), mismatches.other.into_iter()).map(|x| x.flatten()).flatten().collect()
    }
}
