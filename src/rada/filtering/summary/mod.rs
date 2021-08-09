use crate::rada::summarization::{IntervalSummary, LocusSummary};

mod by_editing;

pub use by_editing::SummaryFilterByEditing;

pub trait IntervalSummaryFilter {
    fn is_ok(&self, summary: &IntervalSummary) -> bool;
}

pub trait LocusSummaryFilter {
    fn is_ok(&self, summary: &LocusSummary) -> bool;
}
