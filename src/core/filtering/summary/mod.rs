use crate::core::summary::{LocusSummary, ROISummary};

mod by_mismatches;

pub use by_mismatches::SummaryFilterByMismatches;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait IntervalSummaryFilter {
    fn is_ok(&self, summary: &ROISummary) -> bool;
}

#[cfg_attr(test, automock)]
pub trait LocusSummaryFilter {
    fn is_ok(&self, summary: &LocusSummary) -> bool;
}
