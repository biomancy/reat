use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::ReqNucleotide;
use crate::rada::modules::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};
use crate::rada::modules::summarization::MismatchSummary;

use super::StrandPredictor;

pub struct DummyStrandPredictor {}

impl StrandPredictor for DummyStrandPredictor {}

impl LocusStrandPredictor for DummyStrandPredictor {
    fn predict(&self, _: &Locus, _: &ReqNucleotide, _: &LocusCounts) -> Strand {
        Strand::Unknown
    }
}

impl IntervalStrandPredictor for DummyStrandPredictor {
    fn predict(&self, _: &Interval, _: &MismatchSummary) -> Strand {
        Strand::Unknown
    }
}
