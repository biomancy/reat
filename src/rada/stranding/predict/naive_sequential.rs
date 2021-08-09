// TODO: it should have a static variant based on a macro of some kind
// You create it for 1-12 arguments in a compile time
// AND a from declaration for them using tuples

use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::Nucleotide;
use crate::rada::stranding::predict::{StrandByAtoIEditing, StrandByGenomicFeatures};
use crate::rada::summarization::MismatchesSummary;

use super::{IntervalStrandPredictor, LocusStrandPredictor};

// TODO: Create macro for a generic structure with 1-12 predictors

#[derive(Constructor)]
pub struct NaiveSequentialStrandPredictor {
    by_editing: Option<StrandByAtoIEditing>,
    by_features: Option<StrandByGenomicFeatures>,
}

impl LocusStrandPredictor for NaiveSequentialStrandPredictor {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &LocusCounts) -> Strand {
        let mut strand = Strand::Unknown;

        self.by_features.as_ref().map(|p| strand = LocusStrandPredictor::predict(p, locus, refnuc, sequenced));
        if !strand.is_unknown() {
            return strand;
        }
        self.by_editing.as_ref().map(|p| strand = LocusStrandPredictor::predict(p, locus, refnuc, sequenced));
        strand
    }
}

impl IntervalStrandPredictor for NaiveSequentialStrandPredictor {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand {
        let mut strand = Strand::Unknown;

        self.by_features.as_ref().map(|p| strand = IntervalStrandPredictor::predict(p, interval, mismatches));
        if !strand.is_unknown() {
            return strand;
        }
        self.by_editing.as_ref().map(|p| strand = IntervalStrandPredictor::predict(p, interval, mismatches));
        strand
    }
}
