use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

pub use by_editing::StrandByAtoIEditing;
pub use by_features::StrandByGenomicFeatures;
pub use naive_sequential::NaiveSequentialStrandPredictor;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::Nucleotide;
use crate::rada::summary::MismatchesSummary;

mod by_editing;
mod by_features;
mod naive_sequential;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait IntervalStrandPredictor {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand;
}

impl IntervalStrandPredictor for Box<dyn IntervalStrandPredictor> {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand {
        self.as_ref().predict(interval, mismatches)
    }
}

#[cfg_attr(test, automock)]
pub trait LocusStrandPredictor {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &LocusCounts) -> Strand;
}

impl LocusStrandPredictor for Box<dyn LocusStrandPredictor> {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &LocusCounts) -> Strand {
        self.as_ref().predict(locus, refnuc, sequenced)
    }
}

pub trait StrandPredictor: IntervalStrandPredictor + LocusStrandPredictor {}
