use bio_types::strand::ReqStrand;
use rust_htslib::bam::Record;

mod by_experiment_design;

use crate::rada::read::AlignedRead;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait StrandDeductor<R: AlignedRead> {
    fn deduce(&self, record: &R) -> ReqStrand;
}
