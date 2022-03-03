use bio_types::strand::ReqStrand;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use by_experiment_design::{DeductStrandByDesign, StrandSpecificExperimentDesign};

use crate::core::read::AlignedRead;

mod by_experiment_design;

#[cfg_attr(test, automock)]
pub trait StrandDeductor<R: AlignedRead> {
    fn deduce(&self, record: &R) -> ReqStrand;
}
