// TODO: it should have a static variant based on a macro of some kind
// You create it for 1-12 arguments in a compile time
// AND a from declaration for them using tuples

use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::ReqNucleotide;
use crate::rada::modules::summarization::MismatchSummary;

use super::{IntervalStrandPredictor, LocusStrandPredictor};

pub struct DynamicSequentialStrandPredictor<T> {
    predictors: Vec<T>,
}

impl<T> DynamicSequentialStrandPredictor<T> {
    pub fn new(predictors: Vec<T>) -> Self {
        DynamicSequentialStrandPredictor { predictors }
    }
}

impl<T> LocusStrandPredictor for DynamicSequentialStrandPredictor<T>
where
    T: LocusStrandPredictor,
{
    fn predict(&self, locus: &Locus, refnuc: &ReqNucleotide, sequenced: &LocusCounts) -> Strand {
        for val in &self.predictors {
            let strand = val.predict(locus, refnuc, sequenced);
            if strand != Strand::Unknown {
                return strand;
            }
        }
        Strand::Unknown
    }
}

impl<T> IntervalStrandPredictor for DynamicSequentialStrandPredictor<T>
where
    T: IntervalStrandPredictor,
{
    fn predict(&self, interval: &Interval, mismatches: &MismatchSummary) -> Strand {
        for val in &self.predictors {
            let strand = val.predict(interval, mismatches);
            if strand != Strand::Unknown {
                return strand;
            }
        }
        Strand::Unknown
    }
}
