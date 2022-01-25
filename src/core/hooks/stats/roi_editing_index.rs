use std::fmt::{Display, Formatter};

use bio_types::strand::Strand;
use derive_more::{Add, AddAssign};

use crate::core::hooks::ProcessingHook;
use crate::core::io::Table;
use crate::core::mismatches::roi::{MismatchesSummary, ROIMismatches};

use super::EditingStat;

#[derive(Copy, Clone, Default, Add, AddAssign)]
pub struct ROIEditingIndex {
    accumulator: MismatchesSummary,
}

impl Table for ROIEditingIndex {
    const LENGTH: usize = 16;

    fn row(&self) -> Vec<String> {
        let res = &self.accumulator;
        vec![
            format!("{:.6}", res.A.A as f32 / res.A.coverage() as f32),
            format!("{:.6}", res.T.T as f32 / res.T.coverage() as f32),
            format!("{:.6}", res.G.G as f32 / res.G.coverage() as f32),
            format!("{:.6}", res.C.C as f32 / res.C.coverage() as f32),
            format!("{:.6}", res.A.T as f32 / res.A.coverage() as f32),
            format!("{:.6}", res.T.A as f32 / res.T.coverage() as f32),
            format!("{:.6}", res.A.G as f32 / res.A.coverage() as f32),
            format!("{:.6}", res.T.C as f32 / res.T.coverage() as f32),
            format!("{:.6}", res.A.C as f32 / res.A.coverage() as f32),
            format!("{:.6}", res.T.G as f32 / res.T.coverage() as f32),
            format!("{:.6}", res.G.C as f32 / res.G.coverage() as f32),
            format!("{:.6}", res.C.G as f32 / res.C.coverage() as f32),
            format!("{:.6}", res.G.A as f32 / res.G.coverage() as f32),
            format!("{:.6}", res.C.T as f32 / res.C.coverage() as f32),
            format!("{:.6}", res.G.T as f32 / res.G.coverage() as f32),
            format!("{:.6}", res.C.A as f32 / res.C.coverage() as f32),
        ]
    }

    fn header() -> Vec<String> {
        vec![
            "A->A", "T->T", "G->G", "C->C", "A->T", "T->A", "A->G", "T->C", "A->C", "T->G", "G->C", "C->G", "G->A",
            "C->T", "G->T", "C->A",
        ]
        .into_iter()
        .map(|x| x.to_owned())
        .collect()
    }
}

impl<T: ROIMismatches> ProcessingHook<T> for ROIEditingIndex {
    fn hook(&mut self, objects: Vec<T>) -> Vec<T> {
        for roi in &objects {
            match roi.strand() {
                Strand::Forward => {
                    self.accumulator += *roi.mismatches();
                }
                Strand::Reverse => {
                    self.accumulator += roi.mismatches().complementary();
                }
                Strand::Unknown => {}
            }
        }

        objects
    }
}

impl<T: ROIMismatches> EditingStat<T> for ROIEditingIndex {
    fn combine(stats: &[Self]) -> Self {
        stats.iter().fold(Self::default(), |a, b| a + *b)
    }
}

#[cfg(test)]
mod test {
    use bio_types::genome::Interval;

    use crate::core::dna::NucCounts;
    use crate::core::mismatches::roi::MockROIMismatches;
    use crate::core::workload::ROI;

    use super::*;

    #[test]
    fn header() {
        assert_eq!(
            <ROIEditingIndex as Table>::header(),
            vec![
                "A->A", "T->T", "G->G", "C->C", "A->T", "T->A", "A->G", "T->C", "A->C", "T->G", "G->C", "C->G", "G->A",
                "C->T", "G->T", "C->A"
            ]
        )
    }

    #[test]
    fn row() {
        let factory = |roi, strand, coverage, sequence, mismatches| {
            let mut dummy = MockROIMismatches::new();
            dummy.expect_roi().return_const(roi);
            dummy.expect_strand().return_const(strand);
            dummy.expect_coverage().return_const(coverage);
            dummy.expect_sequence().return_const(sequence);
            dummy.expect_mismatches().return_const(mismatches);
            dummy
        };

        let mut dummy = ROIEditingIndex::default();

        let mut forward = MismatchesSummary::zeros();
        // coverage - 19
        forward.A.A = 10;
        forward.A.C = 2;
        forward.A.G = 3;
        forward.A.T = 4;

        // coverage - 36
        forward.C.A = 15;
        forward.C.C = 6;
        forward.C.G = 7;
        forward.C.T = 8;

        // coverage - 42
        forward.G.A = 9;
        forward.G.C = 10;
        forward.G.G = 11;
        forward.G.T = 12;

        // coverage - 58
        forward.T.A = 13;
        forward.T.C = 14;
        forward.T.G = 15;
        forward.T.T = 16;

        let forward = factory(
            ROI::new(Interval::new("".into(), 1..2), "".into(), vec![1..2]),
            Strand::Forward,
            3u32,
            NucCounts::new(1, 1, 1, 1),
            forward,
        );

        let mut reverse = MismatchesSummary::zeros();
        // coverage - 14
        reverse.A.A = 2;
        reverse.A.C = 3;
        reverse.A.G = 4;
        reverse.A.T = 5;

        // coverage - 30
        reverse.C.A = 6;
        reverse.C.C = 7;
        reverse.C.G = 8;
        reverse.C.T = 9;

        // coverage - 46
        reverse.G.A = 10;
        reverse.G.C = 11;
        reverse.G.G = 12;
        reverse.G.T = 13;

        // coverage - 62
        reverse.T.A = 14;
        reverse.T.C = 15;
        reverse.T.G = 16;
        reverse.T.T = 17;

        let reverse = factory(
            ROI::new(Interval::new("".into(), 1..2), "".into(), vec![1..2]),
            Strand::Reverse,
            9u32,
            NucCounts::new(1, 1, 1, 1),
            reverse,
        );

        let prehook = vec![forward, reverse];
        let posthook = dummy.hook(prehook);

        assert_eq!(
            dummy.row(),
            vec!(
                format!("{:.6}", (10 + 17) as f32 / (19 + 62) as f32),
                format!("{:.6}", (16 + 2) as f32 / (58 + 14) as f32),
                format!("{:.6}", (11 + 7) as f32 / (42 + 30) as f32),
                format!("{:.6}", (6 + 12) as f32 / (36 + 46) as f32),
                format!("{:.6}", (4 + 14) as f32 / (19 + 62) as f32),
                format!("{:.6}", (13 + 5) as f32 / (58 + 14) as f32),
                format!("{:.6}", (3 + 15) as f32 / (19 + 62) as f32),
                format!("{:.6}", (14 + 4) as f32 / (58 + 14) as f32),
                format!("{:.6}", (2 + 16) as f32 / (19 + 62) as f32),
                format!("{:.6}", (15 + 3) as f32 / (58 + 14) as f32),
                format!("{:.6}", (10 + 8) as f32 / (42 + 30) as f32),
                format!("{:.6}", (7 + 11) as f32 / (36 + 46) as f32),
                format!("{:.6}", (9 + 9) as f32 / (42 + 30) as f32),
                format!("{:.6}", (8 + 10) as f32 / (36 + 46) as f32),
                format!("{:.6}", (12 + 6) as f32 / (42 + 30) as f32),
                format!("{:.6}", (15 + 13) as f32 / (36 + 46) as f32),
            )
        )
    }
}
