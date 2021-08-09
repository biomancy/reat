use bio_types::genome::Locus;
use bio_types::strand::Strand;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::ReqNucleotide;

pub struct LocusSummary {
    pub locus: Locus,
    pub strand: Strand,
    pub refnuc: ReqNucleotide,
    pub sequenced: LocusCounts,
}

impl LocusSummary {
    pub fn zeros(locus: Locus, strand: Strand, refnuc: ReqNucleotide) -> Self {
        LocusSummary { locus, strand, refnuc, sequenced: LocusCounts::zeros() }
    }
}
