use std::ops::Range;

use derive_getters::Getters;
use derive_more::Constructor;
use itertools::zip;

use crate::core::hooks::ProcessingHook;
use crate::core::mismatches::interval::{
    IntermediateIntervalMismatches, IntervalMismatches, RefIntervalMismatches, SeparableIntervalMismatches,
};
use crate::core::mismatches::roi::{ROIMismatches, RefROIMismatches};

#[derive(Constructor, Getters, Debug, PartialEq, Copy, Clone)]
pub struct ByMismatches {
    minmismatches: u32,
    minfreq: f32,
    mincov: u32,
}

impl ByMismatches {
    fn ok_roi(&self, x: &impl ROIMismatches) -> bool {
        let (cov, mismatch) = (x.mismatches().coverage(), x.mismatches().mismatches());
        x.coverage() >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }

    fn ok_ranges(&self, x: &impl IntervalMismatches) -> Vec<Range<usize>> {
        let mut result = Vec::new();

        let mut prevind: Option<usize> = None;
        for (ind, (seqnc, refnc)) in zip(x.ncounts(), x.refnuc()).enumerate() {
            let cov = seqnc.coverage();
            let mismatch = seqnc.mismatches(refnc);
            let isok =
                cov >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq;

            match (prevind, isok) {
                // Faulty range continues
                (None, false) => {}
                // Faulty range finished
                (None, true) => prevind = Some(ind),
                // Ok range continues
                (Some(_), true) => {}
                // Ok range finished
                (Some(prev), false) => {
                    result.push(prev..ind);
                    prevind = None;
                }
            }
        }
        if let Some(prev) = prevind {
            result.push(prev..x.ncounts().len());
        }

        result
    }
}

impl<'a> ProcessingHook<RefIntervalMismatches<'a>> for ByMismatches {
    fn hook(&mut self, mut objects: Vec<RefIntervalMismatches<'a>>) -> Vec<RefIntervalMismatches<'a>> {
        let mut result = Vec::new();

        for obj in objects {
            let ranges = self.ok_ranges(&obj);
            obj.extract(&ranges, &mut result);
        }
        result
    }
}

impl<'a> ProcessingHook<RefROIMismatches<'a>> for ByMismatches {
    fn hook(&mut self, mut objects: Vec<RefROIMismatches<'a>>) -> Vec<RefROIMismatches<'a>> {
        objects.retain(|x| self.ok_roi(x));
        objects
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use bio_types::strand::Strand;

    use crate::core::dna::{NucCounts, Nucleotide};
    use crate::core::mismatches::roi::{MismatchesSummary, OwnedROIMismatches};
    use crate::core::workload::ROI;

    use super::*;

    #[test]
    fn ok_roi() {
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
            ROI::new(Interval::new("".into(), 1..2), "".into(), vec![1..2]),
            Strand::Forward,
            0,
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
            assert_eq!(filter.ok_roi(&dummy), expected);
        }
    }

    #[test]
    fn ok_interval() {
        let refnuc = [
            Nucleotide::A,
            Nucleotide::T,
            Nucleotide::G,
            Nucleotide::C,
            Nucleotide::A,
            Nucleotide::Unknown,
            Nucleotide::Unknown,
        ];
        let sequenced = [NucCounts { A: 1, C: 2, G: 3, T: 4 }].repeat(7);

        let dummy = RefIntervalMismatches::new(Interval::new("1".into(), 10..17), Strand::Forward, &refnuc, &sequenced);
        for (expected, minmismatches, minfreq, mincov) in [
            (vec![0..7], 0, 0f32, 0),
            (vec![0..1, 3..7], 8, 0f32, 0),
            (vec![0..1, 4..7], 9, 0.9f32, 9),
            (vec![5..7], 10, 0.95f32, 10),
        ] {
            let filter = ByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(filter.ok_ranges(&dummy), expected);
        }
    }
}
