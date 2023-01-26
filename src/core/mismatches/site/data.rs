use bio_types::genome::Position;
use soa_derive::StructOfArray;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::refpred::PredNucleotide;

#[derive(Clone, Debug, Default, StructOfArray)]
#[soa_derive(Clone, Debug)]
pub struct SiteData {
    // Genomic coordinate for each site
    pub pos: Position,
    // FASTA reference nucleotides
    pub refnuc: Nucleotide,
    // Auto corrected reference nucleotide
    pub prednuc: PredNucleotide,
    // Sequenced nucleotides
    pub sequenced: NucCounts,
}

impl From<SiteDataRef<'_>> for SiteData {
    fn from(x: SiteDataRef<'_>) -> Self {
        Self { pos: *x.pos, refnuc: *x.refnuc, prednuc: *x.prednuc, sequenced: *x.sequenced }
    }
}
