use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::counting::LocusCounts;
use crate::core::dna::Nucleotide;
use crate::core::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};
use crate::core::summary::MismatchesSummary;

use super::StrandPredictor;
use derive_getters::Getters;

#[derive(Constructor, Getters, Copy, Clone)]
pub struct StrandByAtoIEditing {
    min_mismatches: u32,
    min_freq: f32,
}

impl StrandByAtoIEditing {
    #[inline]
    fn is_ok(&self, matches: u32, mismatches: u32) -> bool {
        let coverage = mismatches + matches;
        coverage > 0 && mismatches >= self.min_mismatches && (mismatches as f32 / coverage as f32) >= self.min_freq
    }
}

impl StrandPredictor for StrandByAtoIEditing {}

impl IntervalStrandPredictor for StrandByAtoIEditing {
    fn predict(&self, _: &Interval, mismatches: &MismatchesSummary) -> Strand {
        let a2g = self.is_ok(mismatches.A.A, mismatches.A.G);
        let t2c = self.is_ok(mismatches.T.T, mismatches.T.C);

        match (a2g, t2c) {
            (false, false) => Strand::Unknown,
            (false, true) => Strand::Reverse,
            (true, false) => Strand::Forward,
            (true, true) => {
                let a2g_coverage = mismatches.A.A + mismatches.A.G;
                let t2c_coverage = mismatches.T.T + mismatches.T.C;

                if a2g_coverage == 0 && t2c_coverage == 0 {
                    return Strand::Unknown;
                }

                let a2g = mismatches.A.G as f32 / a2g_coverage as f32;
                let t2c = mismatches.T.C as f32 / t2c_coverage as f32;
                if mismatches.A.G > 0 && a2g > t2c {
                    Strand::Forward
                } else if mismatches.T.C > 0 && t2c > a2g {
                    Strand::Reverse
                } else {
                    Strand::Unknown
                }
            }
        }
    }
}

impl LocusStrandPredictor for StrandByAtoIEditing {
    fn predict(&self, _: &Locus, refnuc: &Nucleotide, sequenced: &LocusCounts) -> Strand {
        match refnuc {
            Nucleotide::A => {
                if self.is_ok(sequenced.A, sequenced.G) {
                    Strand::Forward
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::T => {
                if self.is_ok(sequenced.T, sequenced.C) {
                    Strand::Reverse
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::C | Nucleotide::G | Nucleotide::Unknown => Strand::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use bio_types::strand::Same;

    use super::*;

    #[test]
    fn locuspred() {
        let locus = Locus::new(String::from(""), 123);

        // Not relevant nucleotides
        let a2g = LocusCounts::new(10, 0, 10, 0);
        let t2c = LocusCounts::new(0, 10, 0, 10);
        let dummy = StrandByAtoIEditing::new(1, 0.1);
        for refnuc in [Nucleotide::C, Nucleotide::G, Nucleotide::Unknown] {
            assert!(LocusStrandPredictor::predict(&dummy, &locus, &refnuc, &a2g).is_unknown());
            assert!(LocusStrandPredictor::predict(&dummy, &locus, &refnuc, &t2c).is_unknown());
        }

        // Relevant nucleotides
        let dummy = StrandByAtoIEditing::new(8, 0.05);
        for (result, matches, mismatches) in
            [(&Strand::Forward, 8, 8), (&Strand::Unknown, 100, 4), (&Strand::Unknown, 1, 7), (&Strand::Forward, 10, 10)]
        {
            let a2g = LocusCounts::new(matches, 0, mismatches, 0);
            assert!(LocusStrandPredictor::predict(&dummy, &locus, &Nucleotide::A, &a2g).same(result));

            let t2c = LocusCounts::new(0, mismatches, 0, matches);
            assert!(LocusStrandPredictor::predict(&dummy, &locus, &Nucleotide::T, &t2c).same(&result.neg()));
        }
    }

    #[test]
    fn intervalpred() {
        let interval = Interval::new("".into(), 1..2);
        let dummy = StrandByAtoIEditing::new(8, 0.05);

        // Simple cases
        for (result, matches, mismatches) in
            [(&Strand::Forward, 8, 8), (&Strand::Unknown, 100, 4), (&Strand::Unknown, 1, 7), (&Strand::Forward, 10, 10)]
        {
            let mut summary = MismatchesSummary::zeros();
            summary.A.A = matches;
            summary.A.G = mismatches;
            assert!(IntervalStrandPredictor::predict(&dummy, &interval, &summary).same(&result));

            let mut summary = MismatchesSummary::zeros();
            summary.T.T = matches;
            summary.T.C = mismatches;
            assert!(IntervalStrandPredictor::predict(&dummy, &interval, &summary).same(&result.neg()));
        }

        // Both strands pass the threshold
        for (result, matches, a2g, t2c) in
            [(&Strand::Unknown, 10, 10, 10), (&Strand::Reverse, 10, 10, 11), (&Strand::Forward, 10, 11, 10)]
        {
            let mut summary = MismatchesSummary::zeros();
            summary.A.A = matches;
            summary.A.G = a2g;

            summary.T.T = matches;
            summary.T.C = t2c;

            assert!(IntervalStrandPredictor::predict(&dummy, &interval, &summary).same(&result));
        }
    }
}
