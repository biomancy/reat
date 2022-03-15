#[cfg(test)]
use mockall::{automock, predicate::*};

pub use by_flags::ByFlags;
pub use by_quality::ByQuality;
pub use sequential::Sequential;

use crate::core::read::AlignedRead;

mod by_flags;
mod by_quality;
mod sequential;

#[cfg_attr(test, automock)]
pub trait ReadsFilter<R: AlignedRead> {
    fn is_read_ok(&self, _record: &R) -> bool {
        true
    }
    fn is_base_ok(&self, _record: &R, _base: usize) -> bool {
        true
    }
}
