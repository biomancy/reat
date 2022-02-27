use bio_types::genome::{AbstractLocus, Position};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};

pub use super::SiteMismatches;

#[derive(Clone, Debug, Constructor)]
pub struct REATSiteMismatches {
    contig: String,
    pos: Position,
    strand: Strand,
    refnuc: Nucleotide,
    prednuc: Nucleotide,
    sequenced: NucCounts,
}

impl AbstractLocus for REATSiteMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn pos(&self) -> Position {
        self.pos
    }
}

impl SiteMismatches for REATSiteMismatches {
    fn strand(&self) -> Strand {
        self.strand
    }

    fn refnuc(&self) -> Nucleotide {
        self.refnuc
    }

    fn prednuc(&self) -> Nucleotide {
        self.prednuc
    }

    fn sequenced(&self) -> NucCounts {
        self.sequenced
    }
}
