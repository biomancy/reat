use bio_types::genome::Locus;
use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

#[cfg_attr(test, automock)]
pub trait SiteMismatches {
    fn locus(&self) -> &Locus;
    fn strand(&self) -> Strand;
    fn refnuc(&self) -> Nucleotide;
    fn ncounts(&self) -> &NucCounts;
}
