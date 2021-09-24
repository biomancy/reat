use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use by_editing::StrandByAtoIEditing;
pub use by_features::StrandByGenomicFeatures;
pub use sequential::SequentialStrandPredictor;

use crate::core::summary::{LocusSummary, ROISummary};
use bio_types::genome::{AbstractInterval, AbstractLocus};
use itertools::Itertools;

mod by_editing;
mod by_features;
mod sequential;

#[cfg_attr(test, automock)]
pub trait ROIStrandPredictor {
    fn predict(&self, data: &ROISummary) -> Strand;
    fn batch_predict(&self, data: &[ROISummary]) -> Vec<Strand> {
        debug_assert!(data.iter().map(|x| x.interval.contig()).all_equal());
        debug_assert!(data.windows(2).all(|w| w[0].interval.range().start <= w[1].interval.range().start));
        data.iter().map(|x| self.predict(x)).collect_vec()
    }
}

#[cfg_attr(test, automock)]
pub trait LocusStrandPredictor {
    fn predict(&self, data: &LocusSummary) -> Strand;
    fn batch_predict(&self, data: &[LocusSummary]) -> Vec<Strand> {
        debug_assert!(data.iter().map(|x| x.locus.contig()).all_equal());
        debug_assert!(data.windows(2).all(|w| w[0].locus.pos() <= w[1].locus.pos()));
        data.iter().map(|x| self.predict(x)).collect_vec()
    }
}

pub trait StrandPredictor: ROIStrandPredictor + LocusStrandPredictor {}
