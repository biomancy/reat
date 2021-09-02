use derive_getters::Getters;
use derive_more::Constructor;

use crate::core::dna::Nucleotide;
use crate::core::summary::{LocusSummary, ROISummary};

use super::{IntervalSummaryFilter, LocusSummaryFilter};

#[derive(Constructor, Getters, Debug, PartialEq, Copy, Clone)]
pub struct SummaryFilterByMismatches {
    minmismatches: u32,
    minfreq: f32,
}

impl IntervalSummaryFilter for SummaryFilterByMismatches {
    #[inline]
    fn is_ok(&self, summary: &ROISummary) -> bool {
        let (cov, mismatch) = (summary.mismatches.coverage(), summary.mismatches.mismatches());
        mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }
}

impl LocusSummaryFilter for SummaryFilterByMismatches {
    #[inline]
    fn is_ok(&self, summary: &LocusSummary) -> bool {
        if summary.refnuc == Nucleotide::Unknown {
            // Always pass N's if coverage is sufficient to let the user decide what to do with them
            summary.sequenced.coverage() >= self.minmismatches
        } else {
            let (cov, mismatch) = (summary.sequenced.coverage(), summary.sequenced.mismatches(&summary.refnuc.into()));
            mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
        }
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::{Interval, Locus};
    use bio_types::strand::Strand;

    use crate::core::counting::NucCounts;

    use super::*;

    #[test]
    fn is_ok_interval() {
        let mut dummy = ROISummary {
            interval: Interval::new("".into(), 1..2),
            strand: Strand::Forward,
            name: "".to_string(),
            sequenced: NucCounts::zeros(),
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

        for (expected, minmismatches, minfreq) in
            [(false, 14, 0f32), (true, 13, 0f32), (true, 12, 0f32), (true, 13, 0.48f32), (false, 13, 0.5f32)]
        {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq);
            assert_eq!(IntervalSummaryFilter::is_ok(&filter, &dummy), expected);
        }
    }

    #[test]
    fn is_ok_locus() {
        let mut dummy = LocusSummary::new(
            Locus::new("".into(), 1),
            Strand::Unknown,
            Nucleotide::A,
            NucCounts { A: 1, C: 2, G: 3, T: 4 },
        );

        for (expected, minmismatches, minfreq) in
            [(false, 10, 0f32), (true, 9, 0f32), (true, 8, 0f32), (true, 9, 0.85f32), (false, 9, 0.95f32)]
        {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq);
            assert_eq!(LocusSummaryFilter::is_ok(&filter, &dummy), expected);
        }

        dummy.refnuc = Nucleotide::Unknown;
        for (expected, minmismatches, minfreq) in
            [(true, 10, 0f32), (true, 9, 0f32), (false, 11, 0f32), (true, 10, 1f32), (false, 11, 1f32)]
        {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq);
            assert_eq!(LocusSummaryFilter::is_ok(&filter, &dummy), expected);
        }
    }
}
