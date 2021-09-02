use bio_types::genome::Locus;
use bio_types::strand::Strand;

use crate::core::counting::LocusCounts;
use crate::core::dna::Nucleotide;

pub struct LocusSummary {
    pub locus: Locus,
    pub strand: Strand,
    pub refnuc: Nucleotide,
    pub sequenced: LocusCounts,
}

impl LocusSummary {
    pub fn new(locus: Locus, strand: Strand, refnuc: Nucleotide, sequenced: LocusCounts) -> Self {
        LocusSummary { locus, strand, refnuc, sequenced }
    }
}