use crate::rada::modules::summarization::{IntervalSummary, LocusSummary};

mod byediting;

pub use byediting::SummaryFilterByEditing;

pub trait IntervalSummaryFilter {
    fn is_ok(&self, summary: &IntervalSummary) -> bool;
}

pub trait LocusSummaryFilter {
    fn is_ok(&self, summary: &LocusSummary) -> bool;
}
