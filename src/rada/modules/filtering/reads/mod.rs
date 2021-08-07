#[cfg(test)]
use mockall::{automock, predicate::*};
use rust_htslib::bam::Record;

use crate::rada::modules::read::AlignedRead;
pub use byquality::ReadsFilterByQuality;

mod byquality;

#[cfg_attr(test, automock)]
pub trait ReadsFilter<R: AlignedRead> {
    fn is_read_ok(&self, record: &R) -> bool;
    fn is_base_ok(&self, record: &R, base: usize) -> bool;
}
