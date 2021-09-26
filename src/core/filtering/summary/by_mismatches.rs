use derive_getters::Getters;
use derive_more::Constructor;

use crate::core::summary::{LocusSummary, ROISummary};

use super::{LocusSummaryFilter, ROISummaryFilter};

#[derive(Constructor, Getters, Debug, PartialEq, Copy, Clone)]
pub struct SummaryFilterByMismatches {
    minmismatches: u32,
    minfreq: f32,
    mincov: u32,
}

impl ROISummaryFilter for SummaryFilterByMismatches {
    #[inline]
    fn is_ok(&self, summary: &ROISummary) -> bool {
        let (cov, mismatch) = (summary.mismatches.coverage(), summary.mismatches.mismatches());
        summary.coverage >= self.mincov
            && mismatch >= self.minmismatches
            && mismatch as f32 / cov as f32 >= self.minfreq
    }
}

impl LocusSummaryFilter for SummaryFilterByMismatches {
    #[inline]
    fn is_ok(&self, summary: &LocusSummary) -> bool {
        let cov = summary.sequenced.coverage();
        let mismatch = summary.sequenced.mismatches(&summary.refnuc);
        cov >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::{Interval, Locus};
    use bio_types::strand::Strand;

    use crate::core::counting::NucCounts;

    use super::*;
    use crate::core::dna::Nucleotide;

    #[test]
    fn is_ok_interval() {
        let mut dummy = ROISummary {
            interval: Interval::new("".into(), 1..2),
            strand: Strand::Forward,
            name: "".to_string(),
            coverage: 2,
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

        for (expected, minmismatches, minfreq, mincov) in [
            (false, 14, 0f32, 0),
            (true, 13, 0f32, 0),
            (true, 12, 0f32, 0),
            (true, 13, 0.48f32, 0),
            (false, 13, 0.5f32, 0),
            (true, 13, 0.48f32, 1),
            (true, 13, 0.48f32, 2),
            (false, 13, 0.48f32, 3),
            (false, 13, 0.48f32, 4),
        ] {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(ROISummaryFilter::is_ok(&filter, &dummy), expected);
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

        for (expected, minmismatches, minfreq, mincov) in [
            (false, 10, 0f32, 0),
            (true, 9, 0f32, 5),
            (true, 8, 0f32, 8),
            (true, 9, 0.85f32, 9),
            (false, 9, 0.95f32, 10),
            (true, 9, 0.85f32, 10),
            (false, 9, 0.85f32, 11),
        ] {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(LocusSummaryFilter::is_ok(&filter, &dummy), expected);
        }

        dummy.refnuc = Nucleotide::Unknown;
        for (expected, minmismatches, minfreq, mincov) in [
            (true, 10, 0f32, 0),
            (true, 9, 0f32, 0),
            (false, 11, 0f32, 0),
            (true, 10, 1f32, 0),
            (false, 11, 1f32, 0),
            (true, 10, 1f32, 10),
            (false, 10, 1f32, 11),
        ] {
            let filter = SummaryFilterByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(LocusSummaryFilter::is_ok(&filter, &dummy), expected);
        }
    }
}
