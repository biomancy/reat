use bio_types::genome::Locus;
use bio_types::strand::Strand;

use crate::core::counting::NucCounts;
use crate::core::dna::Nucleotide;

pub struct LocusSummary {
    pub locus: Locus,
    pub strand: Strand,
    pub refnuc: Nucleotide,
    pub sequenced: NucCounts,
}

impl LocusSummary {
    pub fn new(locus: Locus, strand: Strand, refnuc: Nucleotide, sequenced: NucCounts) -> Self {
        LocusSummary { locus, strand, refnuc, sequenced }
    }
}
