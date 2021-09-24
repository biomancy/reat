pub use by_heuristic::RefNucPredByHeurisitc;

use crate::core::dna::Nucleotide;

use crate::core::counting::NucCounts;

#[cfg(test)]
use mockall::{automock, predicate::*};

mod by_heuristic;

#[cfg_attr(test, automock)]
pub trait RefNucPredictor {
    fn predict(&self, assembly: &Nucleotide, sequenced: &NucCounts) -> Nucleotide;
}
