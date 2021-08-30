pub use by_heuristic::RefNucPredByHeurisitc;

use crate::core::dna::Nucleotide;

use crate::core::counting::LocusCounts;

#[cfg(test)]
use mockall::{automock, predicate::*};

mod by_heuristic;

#[cfg_attr(test, automock)]
pub trait RefNucPredictor {
    fn predict(&self, assembly: &Nucleotide, sequenced: &LocusCounts) -> Nucleotide;
}
