pub use base::BaseSummarizer;

use crate::rada::counting::NucCounterContent;
use crate::rada::summarization::summary::{IntervalSummary, LocusSummary};

mod base;

pub trait IntervalSummarizer {
    fn summarize(&self, name: &str, content: NucCounterContent) -> Vec<IntervalSummary>;
}

pub trait LocusSummarizer {
    fn summarize(&self, content: NucCounterContent) -> Vec<LocusSummary>;
}

pub trait Summarizer: IntervalSummarizer + LocusSummarizer {}
