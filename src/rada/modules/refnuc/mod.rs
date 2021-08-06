pub use byheuristic::RefNucPredByHeurisitc;

use crate::rada::modules::dna::{Nucleotide, ReqNucleotide};

use super::counting::LocusCounts;

mod byheuristic;

pub trait RefNucPredictor {
    // Only for the + strand! Counts also for the + strand!
    fn predict(
        &self,
        assembly: Nucleotide,
        variants: &[ReqNucleotide],
        sequenced: &LocusCounts,
    ) -> ReqNucleotide;
}
