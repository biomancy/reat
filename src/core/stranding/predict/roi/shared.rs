use bio_types::strand::Strand;

use crate::core::mismatches::roi::IntermediateROIMismatches;
use crate::core::stranding::predict::shared::StrandByAtoIEditing;
use crate::core::stranding::predict::shared::StrandByGenomicAnnotation;

use super::ROIStrandPredictor;

impl<T: IntermediateROIMismatches> ROIStrandPredictor<T> for StrandByAtoIEditing {
    fn predict(&self, mut rois: Vec<T>) -> Vec<T> {
        for x in &mut rois {
            let mismatches = x.mismatches();
            let a2g = self.edited(mismatches.A.A, mismatches.A.G);
            let t2c = self.edited(mismatches.T.T, mismatches.T.C);

            let strand = match (a2g, t2c) {
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
            };
            x.set_strand(strand);
        }
        rois
    }
}

impl<T: IntermediateROIMismatches> ROIStrandPredictor<T> for StrandByGenomicAnnotation {
    fn predict(&self, mut rois: Vec<T>) -> Vec<T> {
        for r in &mut rois {
            r.set_strand(self.predict(r.interval()));
        }
        rois
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use bio_types::genome::Interval;
    use bio_types::genome::Locus;
    use bio_types::strand::Same;
    use itertools::{zip, Itertools};
    use mockall::predicate::function;

    use crate::core::dna::NucCounts;
    use crate::core::mismatches::roi::{MismatchesSummary, MockROIMismatches};

    use super::*;

    #[test]
    fn roi_strand_by_editing() {
        let dummy = StrandByAtoIEditing::new(8, 0.05);

        let mockit = |mismatches, expected| {
            let mut mock = MockROIMismatches::new();
            mock.expect_mismatches().return_const(mismatches);
            mock.expect_set_strand().once().with(function(move |x: &Strand| x.same(&expected))).return_const(());
            mock
        };

        let mut mocks = Vec::new();
        for (result, matches, mismatches) in
            [(Strand::Forward, 8, 8), (Strand::Unknown, 100, 4), (Strand::Unknown, 1, 7), (Strand::Forward, 10, 10)]
        {
            let mut summary = MismatchesSummary::zeros();
            summary.T.T = matches;
            summary.T.C = mismatches;
            mocks.push(mockit(summary, result.clone().neg()));

            let mut summary = MismatchesSummary::zeros();
            summary.A.A = matches;
            summary.A.G = mismatches;
            mocks.push(mockit(summary, result));
        }
        let mut mocks = ROIStrandPredictor::predict(&dummy, mocks);
        mocks.iter_mut().for_each(|x| x.checkpoint());
        mocks.clear();

        for (result, matches, a2g, t2c) in
            [(Strand::Unknown, 10, 10, 10), (Strand::Reverse, 10, 10, 11), (Strand::Forward, 10, 11, 10)]
        {
            let mut summary = MismatchesSummary::zeros();
            summary.A.A = matches;
            summary.A.G = a2g;

            summary.T.T = matches;
            summary.T.C = t2c;

            mocks.push(mockit(summary, result));
        }
        let mut mocks = ROIStrandPredictor::predict(&dummy, mocks);
        mocks.iter_mut().for_each(|x| x.checkpoint());
    }
}
