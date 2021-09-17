use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

pub use by_editing::StrandByAtoIEditing;
pub use by_features::StrandByGenomicFeatures;
pub use sequential::SequentialStrandPredictor;

use crate::core::counting::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::summary::MismatchesSummary;

mod by_editing;
mod by_features;
mod sequential;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait ROIStrandPredictor {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand;
}

impl ROIStrandPredictor for Box<dyn ROIStrandPredictor> {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand {
        self.as_ref().predict(interval, mismatches)
    }
}

#[cfg_attr(test, automock)]
pub trait LocusStrandPredictor {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &NucCounts) -> Strand;
}

impl LocusStrandPredictor for Box<dyn LocusStrandPredictor> {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &NucCounts) -> Strand {
        self.as_ref().predict(locus, refnuc, sequenced)
    }
}

pub trait StrandPredictor: ROIStrandPredictor + LocusStrandPredictor {}
