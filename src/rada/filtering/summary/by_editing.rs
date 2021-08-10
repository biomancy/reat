use derive_more::Constructor;

use crate::rada::summary::IntervalSummary;

use super::IntervalSummaryFilter;

#[derive(Constructor, Debug, PartialEq, Clone)]
pub struct SummaryFilterByEditing {
    minmismatches: u32,
    minfreq: f32,
}

impl IntervalSummaryFilter for SummaryFilterByEditing {
    fn is_ok(&self, summary: &IntervalSummary) -> bool {
        let (cov, mismatch) = (summary.mismatches.coverage(), summary.mismatches.mismatches());
        mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use bio_types::strand::Strand;

    use super::*;

    #[test]
    fn is_ok_by_mismatches() {
        let mut dummy = IntervalSummary {
            interval: Interval::new("".to_string(), 1..2),
            strand: Strand::Forward,
            name: "".to_string(),
            mismatches: Default::default(),
        };

        dummy.mismatches.A.C = 1;
        dummy.mismatches.A.A = 4;

        dummy.mismatches.C.G = 1;
        dummy.mismatches.C.C = 2;

        dummy.mismatches.G.T = 1;
        dummy.mismatches.G.G = 5;

        dummy.mismatches.T.A = 10;
        dummy.mismatches.T.T = 3;

        assert!(!SummaryFilterByEditing::new(14, 0f32).is_ok(&dummy));
        assert!(SummaryFilterByEditing::new(13, 0f32).is_ok(&dummy));
        assert!(SummaryFilterByEditing::new(12, 0f32).is_ok(&dummy));
    }
}
