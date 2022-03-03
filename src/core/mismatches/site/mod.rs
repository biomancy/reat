use std::ops::Range;

use bio_types::genome::{AbstractLocus, Position};
use bio_types::strand::Strand;

pub use batched::REATBatchedSiteMismatches;
pub use builder::REATSiteMismatchesBuilder;
pub use flat::REATSiteMismatches;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

use super::BatchedMismatches;

mod batched;
mod builder;
mod flat;

pub type SiteMismatchesPreview = (Nucleotide, NucCounts);

pub trait SiteMismatches: AbstractLocus {
    fn trstrand(&self) -> Strand;
    fn refnuc(&self) -> Nucleotide;
    fn prednuc(&self) -> Nucleotide;
    fn sequenced(&self) -> NucCounts;
}

pub trait BatchedSiteMismatches: BatchedMismatches
    where
        Self::Flattened: SiteMismatches,
{
    fn pos(&self) -> &[Position];
    fn refnuc(&self) -> &[Nucleotide];
    fn prednuc(&self) -> &[Nucleotide];
    fn seqnuc(&self) -> &[NucCounts];
    // Continuous blocks in the genome coordinates
    fn blocks(&self) -> &[Range<Position>];
}
