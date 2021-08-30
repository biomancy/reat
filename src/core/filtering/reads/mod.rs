#[cfg(test)]
use mockall::{automock, predicate::*};

pub use by_quality::ReadsFilterByQuality;

use crate::core::read::AlignedRead;

mod by_quality;

#[cfg_attr(test, automock)]
pub trait ReadsFilter<R: AlignedRead> {
    fn is_read_ok(&self, record: &R) -> bool;
    fn is_base_ok(&self, record: &R, base: usize) -> bool;
}
