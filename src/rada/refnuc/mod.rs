pub use by_heuristic::RefNucPredByHeurisitc;

use crate::rada::dna::Nucleotide;

use crate::rada::counting::LocusCounts;

#[cfg(test)]
use mockall::{automock, predicate::*};

mod by_heuristic;

#[cfg_attr(test, automock)]
pub trait RefNucPredictor {
    fn predict(&self, assembly: &Nucleotide, sequenced: &LocusCounts) -> Nucleotide;
}
