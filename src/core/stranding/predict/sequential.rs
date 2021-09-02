// TODO: it should have a static variant based on a macro of some kind
// You create it for 1-12 arguments in a compile time
// AND a from declaration for them using tuples

use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::counting::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::stranding::predict::{StrandByAtoIEditing, StrandByGenomicFeatures};
use crate::core::summary::MismatchesSummary;

use super::{IntervalStrandPredictor, LocusStrandPredictor};

#[derive(Constructor, Clone)]
pub struct SequentialStrandPredictor {
    by_editing: Option<StrandByAtoIEditing>,
    by_features: Option<StrandByGenomicFeatures>,
}

impl SequentialStrandPredictor {
    pub fn by_editing(&self) -> &Option<StrandByAtoIEditing> {
        &self.by_editing
    }
    pub fn by_features(&self) -> &Option<StrandByGenomicFeatures> {
        &self.by_features
    }
}

impl LocusStrandPredictor for SequentialStrandPredictor {
    fn predict(&self, locus: &Locus, refnuc: &Nucleotide, sequenced: &NucCounts) -> Strand {
        let mut strand = Strand::Unknown;

        if let Some(p) = self.by_features.as_ref() {
            strand = LocusStrandPredictor::predict(p, locus, refnuc, sequenced);
        }
        if !strand.is_unknown() {
            return strand;
        }
        if let Some(p) = self.by_editing.as_ref() {
            strand = LocusStrandPredictor::predict(p, locus, refnuc, sequenced);
        }
        strand
    }
}

impl IntervalStrandPredictor for SequentialStrandPredictor {
    fn predict(&self, interval: &Interval, mismatches: &MismatchesSummary) -> Strand {
        let mut strand = Strand::Unknown;

        if let Some(p) = self.by_features.as_ref() {
            strand = IntervalStrandPredictor::predict(p, interval, mismatches);
        }
        if !strand.is_unknown() {
            return strand;
        }
        if let Some(p) = self.by_editing.as_ref() {
            strand = IntervalStrandPredictor::predict(p, interval, mismatches);
        }
        strand
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;
//
//     #[test]
//     fn predict() {
//         // TODO: make a reasonable unit test once Sequential predictor is a macro
//     }
// }
