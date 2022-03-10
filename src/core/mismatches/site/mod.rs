use std::ops::Range;

use bio_types::genome::{AbstractLocus, Position};
use bio_types::strand::Strand;

pub use builder::REATSiteMismatchesBuilder;
pub use flat::REATSiteMismatches;
pub use vec::REATSiteMismatchesVec;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

use super::MismatchesVec;

mod builder;
mod flat;
mod vec;

pub type SiteMismatchesPreview = (Nucleotide, NucCounts);

pub trait SiteMismatches: AbstractLocus {
    // Transcription strand
    fn trstrand(&self) -> Strand;
    // Reference nucleotide (FASTA)
    fn refnuc(&self) -> Nucleotide;
    // Autoref-corrected reference
    fn prednuc(&self) -> Nucleotide;
    // Sequenced nucleotides
    fn sequenced(&self) -> NucCounts;
}

pub trait SiteMismatchesVec: MismatchesVec
where
    Self::Flat: SiteMismatches,
{
    // Genomic coordinate for each loci
    fn pos(&self) -> &[Position];
    // FASTA reference nucleotides
    fn refnuc(&self) -> &[Nucleotide];
    // Autoref-corrected reference nucleotide
    fn prednuc(&self) -> &[Nucleotide];
    // Sequenced nucleotides
    fn seqnuc(&self) -> &[NucCounts];
    // Continuous blocks in the genome coordinates
    fn blocks(&self) -> &[Range<Position>];
}
