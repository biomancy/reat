use bio_types::genome::Locus;
use bio_types::strand::Strand;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::{Nucleotide, ReqNucleotide};

pub struct LocusSummary {
    pub locus: Locus,
    pub strand: Strand,
    pub refnuc: Nucleotide,
    pub sequenced: LocusCounts,
}

impl LocusSummary {
    pub fn zeros(locus: Locus, strand: Strand, refnuc: Nucleotide) -> Self {
        LocusSummary { locus, strand, refnuc, sequenced: LocusCounts::zeros() }
    }

    pub fn new(locus: Locus, strand: Strand, refnuc: Nucleotide, sequenced: LocusCounts) -> Self {
        LocusSummary { locus, strand, refnuc, sequenced }
    }
}
