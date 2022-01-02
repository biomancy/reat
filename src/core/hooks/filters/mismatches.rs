use derive_getters::Getters;
use derive_more::Constructor;
use itertools::zip;

use crate::core::hooks::ProcessingHook;
use crate::core::mismatches::interval::{BorrowedIntervalMismatches, IntervalMismatches};
use crate::core::mismatches::roi::{BorrowedROIMismatches, ROIMismatches};

#[derive(Constructor, Getters, Debug, PartialEq, Copy, Clone)]
pub struct ByMismatches {
    minmismatches: u32,
    minfreq: f32,
    mincov: u32,
}

impl ByMismatches {
    fn is_roi_ok(&self, x: &impl ROIMismatches) -> bool {
        let (cov, mismatch) = (x.mismatches().coverage(), x.mismatches().mismatches());
        x.coverage() >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }

    fn is_interval_ok(&self, x: &impl IntervalMismatches) -> Vec<bool> {
        zip(x.ncounts(), x.refnuc())
            .map(|(seqnc, refnc)| {
                let cov = seqnc.coverage();
                let mismatch = seqnc.mismatches(refnc);
                cov >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
            })
            .collect()
    }
}

impl<'a> ProcessingHook<BorrowedIntervalMismatches<'a>> for ByMismatches {
    fn hook(&mut self, objects: Vec<BorrowedIntervalMismatches<'a>>) -> Vec<BorrowedIntervalMismatches<'a>> {
        todo!()
    }
}

impl<'a> ProcessingHook<BorrowedROIMismatches<'a>> for ByMismatches {
    fn hook(&mut self, objects: Vec<BorrowedROIMismatches<'a>>) -> Vec<BorrowedROIMismatches<'a>> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use bio_types::strand::Strand;

    use crate::core::dna::NucCounts;
    use crate::core::mismatches::roi::{MismatchesSummary, OwnedROIMismatches};

    use super::*;

    #[test]
    fn is_roi_ok() {
        let mut dummy: MismatchesSummary = Default::default();

        dummy.A.C = 1;
        dummy.A.A = 4;

        dummy.C.G = 1;
        dummy.C.C = 2;

        dummy.G.T = 1;
        dummy.G.G = 5;

        dummy.T.A = 10;
        dummy.T.T = 3;

        let mut dummy = OwnedROIMismatches::new(
            Interval::new("".into(), 1..2),
            Strand::Forward,
            "".to_string(),
            2,
            NucCounts::zeros(),
            dummy,
        );

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
            let filter = ByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(filter.is_roi_ok(&dummy), expected);
        }
    }

    // #[test]
    // fn is_ok_locus() {
    //     let mut dummy = BlockSiteSummary {
    //         block: Interval::new("".into(), 0..1),
    //         strand: Strand::Forward,
    //         refnuc: vec![Nucleotide::A],
    //         ncounts: vec![NucCounts { A: 1, C: 2, G: 3, T: 4 }],
    //     };
    //
    //     for (expected, minmismatches, minfreq, mincov) in [
    //         (false, 10, 0f32, 0),
    //         (true, 9, 0f32, 5),
    //         (true, 8, 0f32, 8),
    //         (true, 9, 0.85f32, 9),
    //         (false, 9, 0.95f32, 10),
    //         (true, 9, 0.85f32, 10),
    //         (false, 9, 0.85f32, 11),
    //     ] {
    //         let filter = SitesFilteringByMismatches::new(minmismatches, minfreq, mincov);
    //         assert_eq!(filter.isok(&dummy)[0], expected);
    //     }
    //
    //     dummy.refnuc[0] = Nucleotide::Unknown;
    //     for (expected, minmismatches, minfreq, mincov) in [
    //         (true, 10, 0f32, 0),
    //         (true, 9, 0f32, 0),
    //         (false, 11, 0f32, 0),
    //         (true, 10, 1f32, 0),
    //         (false, 11, 1f32, 0),
    //         (true, 10, 1f32, 10),
    //         (false, 10, 1f32, 11),
    //     ] {
    //         let filter = SitesFilteringByMismatches::new(minmismatches, minfreq, mincov);
    //         assert_eq!(filter.isok(&dummy)[0], expected);
    //     }
    // }
}
