use bio_types::strand::{Same, Strand};
use derive_getters::Getters;
use derive_more::Constructor;
use itertools::zip;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{BinnedROIMismatches, MismatchesSummary, REATBinnedROIMismatches};
use crate::core::mismatches::site::{BinnedSiteMismatches, REATBatchedSiteMismatches};
use crate::core::mismatches::BatchedMismatches;
use crate::core::stranding::predict::ctx::StrandingAlgoResult;
use crate::core::stranding::predict::{StrandingAlgo, StrandingContext};

#[derive(Constructor, Getters, Copy, Clone)]
pub struct StrandByAtoIEditing {
    minmismatches: u32,
    minfreq: f32,
}

impl StrandByAtoIEditing {
    #[inline]
    pub fn edited(&self, matches: u32, mismatches: u32) -> bool {
        let coverage = mismatches + matches;
        coverage > 0 && mismatches >= self.minmismatches && (mismatches as f32 / coverage as f32) >= self.minfreq
    }

    #[inline]
    fn nucpred(&self, sequenced: &NucCounts, refnuc: &Nucleotide) -> Strand {
        match refnuc {
            Nucleotide::A => {
                if self.edited(sequenced.A, sequenced.G) {
                    Strand::Forward
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::T => {
                if self.edited(sequenced.T, sequenced.C) {
                    Strand::Reverse
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::C | Nucleotide::G | Nucleotide::Unknown => Strand::Unknown,
        }
    }

    #[inline]
    fn roipred(&self, mismatches: &MismatchesSummary) -> Strand {
        let a2g = self.edited(mismatches.A.A, mismatches.A.G);
        let t2c = self.edited(mismatches.T.T, mismatches.T.C);

        match (a2g, t2c) {
            (false, false) => Strand::Unknown,
            (false, true) => Strand::Reverse,
            (true, false) => Strand::Forward,
            (true, true) => {
                let a2g_coverage = mismatches.A.A + mismatches.A.G;
                let t2c_coverage = mismatches.T.T + mismatches.T.C;

                if a2g_coverage == 0 && t2c_coverage == 0 {
                    Strand::Unknown
                } else {
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
}

impl StrandingAlgo<REATBinnedROIMismatches> for StrandByAtoIEditing {
    fn predict(&self, mut ctx: StrandingContext<REATBinnedROIMismatches>) -> StrandingContext<REATBinnedROIMismatches> {
        ctx.apply(|x| StrandingAlgoResult::EachElement((x.mismatches().iter().map(|m| self.roipred(m)).collect())));
        ctx
    }
}

impl StrandingAlgo<REATBatchedSiteMismatches> for StrandByAtoIEditing {
    fn predict(
        &self,
        mut ctx: StrandingContext<REATBatchedSiteMismatches>,
    ) -> StrandingContext<REATBatchedSiteMismatches> {
        ctx.apply(|x| {
            StrandingAlgoResult::EachElement(
                (zip(x.prednuc(), x.seqnuc()).map(|(pn, sn)| self.nucpred(sn, pn))).collect(),
            )
        });
        ctx
    }
}

// #[cfg(test)]
// mod tests {
//     use std::ops::Neg;
//
//     use bio_types::genome::Locus;
//     use bio_types::genome::{AbstractInterval, Interval};
//     use bio_types::strand::Same;
//     use bio_types::strand::Strand;
//     use itertools::izip;
//     use itertools::{zip, Itertools};
//     use mockall::predicate::function;
//
//     use crate::core::dna::{NucCounts, Nucleotide};
//     use crate::core::mismatches::interval::{IntervalMismatches, MockIntervalMismatches, RefIntervalMismatches};
//     use crate::core::mismatches::roi::{MismatchesSummary, MockROIMismatches};
//
//     use super::*;
//
//     #[test]
//     fn roi_strand_by_editing() {
//         let dummy = StrandByAtoIEditing::new(8, 0.05);
//
//         let mockit = |mismatches, expected| {
//             let mut mock = MockROIMismatches::new();
//             mock.expect_mismatches().return_const(mismatches);
//             mock.expect_set_strand().once().with(function(move |x: &Strand| x.same(&expected))).return_const(());
//             mock
//         };
//
//         let mut mocks = Vec::new();
//         for (result, matches, mismatches) in
//             [(Strand::Forward, 8, 8), (Strand::Unknown, 100, 4), (Strand::Unknown, 1, 7), (Strand::Forward, 10, 10)]
//         {
//             let mut summary = MismatchesSummary::zeros();
//             summary.T.T = matches;
//             summary.T.C = mismatches;
//             mocks.push(mockit(summary, result.clone().neg()));
//
//             let mut summary = MismatchesSummary::zeros();
//             summary.A.A = matches;
//             summary.A.G = mismatches;
//             mocks.push(mockit(summary, result));
//         }
//         let mut mocks = ROIStrandPredictor::predict(&dummy, mocks);
//         mocks.iter_mut().for_each(|x| x.checkpoint());
//         mocks.clear();
//
//         for (result, matches, a2g, t2c) in
//             [(Strand::Unknown, 10, 10, 10), (Strand::Reverse, 10, 10, 11), (Strand::Forward, 10, 11, 10)]
//         {
//             let mut summary = MismatchesSummary::zeros();
//             summary.A.A = matches;
//             summary.A.G = a2g;
//
//             summary.T.T = matches;
//             summary.T.C = t2c;
//
//             mocks.push(mockit(summary, result));
//         }
//         let mut mocks = ROIStrandPredictor::predict(&dummy, mocks);
//         mocks.iter_mut().for_each(|x| x.checkpoint());
//     }
//
//     #[test]
//     fn nucpred() {
//         let dummy = StrandByAtoIEditing::new(10, 0.1);
//
//         // Unknown strand
//         for (sequenced, refnuc) in [
//             (NucCounts::T(123), Nucleotide::C),
//             (NucCounts::G(234), Nucleotide::Unknown),
//             (NucCounts::A(32), Nucleotide::G),
//             (NucCounts::new(170, 0, 170, 0), Nucleotide::C),
//             (NucCounts::new(10, 0, 9, 0), Nucleotide::A),
//             (NucCounts::new(0, 200, 0, 200), Nucleotide::G),
//         ] {
//             assert!(dummy.nucpred(&sequenced, &refnuc).is_unknown());
//         }
//
//         let dummy = StrandByAtoIEditing::new(8, 0.05);
//
//         // Inferred strand
//         for (matches, mismatches, strand) in
//             [(8, 8, Strand::Forward), (100, 4, Strand::Unknown), (1, 7, Strand::Unknown), (10, 10, Strand::Forward)]
//         {
//             let cnts = NucCounts::new(matches, 0, mismatches, 0);
//             assert!(dummy.nucpred(&cnts, &Nucleotide::A).same(&strand));
//
//             let cnts = NucCounts::new(0, mismatches, 0, matches);
//             assert!(dummy.nucpred(&cnts, &Nucleotide::T).same(&strand.neg()));
//         }
//     }
//
//     #[test]
//     fn intervals_by_editing() {
//         let dummy = StrandByAtoIEditing::new(10, 0.5);
//         let interval = Interval::new("".into(), 100..101);
//
//         let assertok = |workload: Vec<(Vec<NucCounts>, Vec<Nucleotide>)>, strands, explen| {
//             let mut blocks = Vec::new();
//             for (sequenced, reference) in &workload {
//                 blocks.push(RefIntervalMismatches::new(interval.clone(), Strand::Unknown, &reference, &sequenced))
//             }
//
//             let predicted = IntervalStrandPredictor::predict(&dummy, blocks);
//             for (pred, expstrnd, explen) in izip!(predicted, strands, explen) {
//                 assert_eq!(pred.ncounts().len(), pred.refnuc().len());
//                 assert_eq!(pred.refnuc().len(), explen);
//                 assert!(pred.strand().same(&expstrnd));
//             }
//         };
//
//         // Empty list
//         let result: Vec<MockIntervalMismatches> = IntervalStrandPredictor::predict(&dummy, vec![]);
//         assert!(result.is_empty());
//
//         // Dummmies
//         let a2g = |len| (vec![NucCounts { A: 5, C: 1, G: 20, T: 3 }].repeat(len), vec![Nucleotide::A].repeat(len));
//         let t2c = |len| (vec![NucCounts { A: 2, C: 235, G: 4, T: 10 }].repeat(len), vec![Nucleotide::T].repeat(len));
//         let unknown = |len| (vec![NucCounts { A: 1, C: 2, G: 3, T: 4 }].repeat(len), vec![Nucleotide::G].repeat(len));
//         let chained = |items: Vec<(Vec<NucCounts>, Vec<Nucleotide>)>| {
//             let mut result = (Vec::new(), Vec::new());
//             for item in items {
//                 result.0.extend(item.0);
//                 result.1.extend(item.1);
//             }
//             result
//         };
//
//         // Single nucleotide
//         assertok(
//             vec![a2g(1), t2c(1), unknown(1)],
//             vec![Strand::Forward, Strand::Reverse, Strand::Unknown],
//             vec![1, 1, 1],
//         );
//
//         // Interval of nucleotides with each possible strand
//         assertok(
//             vec![a2g(10), t2c(3), unknown(150)],
//             vec![Strand::Forward, Strand::Reverse, Strand::Unknown],
//             vec![10, 3, 150],
//         );
//
//         // Intervals of nucleotides with multiple strands
//         assertok(
//             vec![
//                 chained(vec![a2g(13), t2c(49), a2g(1)]),
//                 chained(vec![a2g(1), t2c(12), unknown(199)]),
//                 chained(vec![t2c(123), unknown(1), a2g(1)]),
//             ],
//             vec![
//                 Strand::Forward,
//                 Strand::Reverse,
//                 Strand::Forward,
//                 Strand::Forward,
//                 Strand::Reverse,
//                 Strand::Unknown,
//                 Strand::Reverse,
//                 Strand::Unknown,
//                 Strand::Forward,
//             ],
//             vec![13, 49, 1, 1, 12, 199, 123, 1, 1],
//         )
//     }
// }
