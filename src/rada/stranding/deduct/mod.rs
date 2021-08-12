use bio_types::strand::ReqStrand;

mod by_experiment_design;

use crate::rada::read::AlignedRead;
pub use by_experiment_design::{DeductStrandByDesign, StrandSpecificExperimentDesign};

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait StrandDeductor<R: AlignedRead> {
    fn deduce(&self, record: &R) -> ReqStrand;
}
