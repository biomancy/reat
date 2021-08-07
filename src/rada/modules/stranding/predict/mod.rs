use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

pub use byediting::StrandByAtoIEditing;
pub use byfeatures::StrandByGenomicFeatures;
pub use sequential::DynamicSequentialStrandPredictor;

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::ReqNucleotide;
use crate::rada::modules::summarization::MismatchSummary;

mod byediting;
mod byfeatures;
mod dummy;
mod sequential;

pub trait IntervalStrandPredictor {
    fn predict(&self, interval: &Interval, mismatches: &MismatchSummary) -> Strand;
}

impl IntervalStrandPredictor for Box<dyn IntervalStrandPredictor> {
    fn predict(&self, interval: &Interval, mismatches: &MismatchSummary) -> Strand {
        self.as_ref().predict(interval, mismatches)
    }
}

pub trait LocusStrandPredictor {
    fn predict(&self, locus: &Locus, refnuc: &ReqNucleotide, sequenced: &LocusCounts) -> Strand;
}

impl LocusStrandPredictor for Box<dyn LocusStrandPredictor> {
    fn predict(&self, locus: &Locus, refnuc: &ReqNucleotide, sequenced: &LocusCounts) -> Strand {
        self.as_ref().predict(locus, refnuc, sequenced)
    }
}

pub trait StrandPredictor: IntervalStrandPredictor + LocusStrandPredictor {}
