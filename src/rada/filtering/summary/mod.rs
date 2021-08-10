use crate::rada::summary::{IntervalSummary, LocusSummary};

mod by_editing;

pub use by_editing::SummaryFilterByEditing;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait IntervalSummaryFilter {
    fn is_ok(&self, summary: &IntervalSummary) -> bool;
}

#[cfg_attr(test, automock)]
pub trait LocusSummaryFilter {
    fn is_ok(&self, summary: &LocusSummary) -> bool;
}
