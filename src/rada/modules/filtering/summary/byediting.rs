use crate::rada::modules::summarization::IntervalSummary;

use super::{IntervalSummaryFilter, LocusSummaryFilter};
use derive_more::Constructor;

#[derive(Constructor, Debug, PartialEq, Clone)]
pub struct SummaryFilterByEditing {
    minmismatches: u32,
    minfreq: f32,
}

impl IntervalSummaryFilter for SummaryFilterByEditing {
    fn is_ok(&self, summary: &IntervalSummary) -> bool {
        let (cov, mismatch) = (
            summary.mismatches.coverage(),
            summary.mismatches.mismatches(),
        );
        mismatch > self.minmismatches && mismatch as f32 / cov as f32 > self.minfreq
    }
}
