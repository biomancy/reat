// TODO: it should have a static variant based on a macro of some kind
// You create it for 1-12 arguments in a compile time
// AND a from declaration for them using tuples

use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::stranding::predict::{StrandByAtoIEditing, StrandByGenomicFeatures};
use crate::core::summary::{LocusSummary, ROISummary};

use super::{LocusStrandPredictor, ROIStrandPredictor};
use itertools::zip;

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
    fn predict(&self, data: &LocusSummary) -> Strand {
        let mut strand = Strand::Unknown;

        if let Some(p) = self.by_features.as_ref() {
            strand = LocusStrandPredictor::predict(p, data);
        }
        if !strand.is_unknown() {
            return strand;
        }
        if let Some(p) = self.by_editing.as_ref() {
            strand = LocusStrandPredictor::predict(p, data);
        }
        strand
    }

    fn batch_predict(&self, data: &[LocusSummary]) -> Vec<Strand> {
        if data.is_empty() {
            return vec![];
        }

        let mut strands = if let Some(p) = self.by_features.as_ref() {
            LocusStrandPredictor::batch_predict(p, data)
        } else {
            let mut v = Vec::with_capacity(data.len());
            v.resize(data.len(), Strand::Unknown);
            v
        };

        if let Some(p) = self.by_editing.as_ref() {
            for (s, d) in zip(&mut strands, data) {
                if s.is_unknown() {
                    *s = LocusStrandPredictor::predict(p, d);
                }
            }
        }

        strands
    }
}

impl ROIStrandPredictor for SequentialStrandPredictor {
    fn predict(&self, data: &ROISummary) -> Strand {
        let mut strand = Strand::Unknown;

        if let Some(p) = self.by_features.as_ref() {
            strand = ROIStrandPredictor::predict(p, data);
        }
        if !strand.is_unknown() {
            return strand;
        }
        if let Some(p) = self.by_editing.as_ref() {
            strand = ROIStrandPredictor::predict(p, data);
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
