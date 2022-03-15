use bio_types::strand::ReqStrand;

pub use by_experiment_design::{DeduceStrandByDesign, StrandSpecificExperimentDesign};

use crate::core::read::AlignedRead;

mod by_experiment_design;

pub trait StrandDeducer<R: AlignedRead> {
    fn deduce(&self, record: &R) -> ReqStrand;
}
