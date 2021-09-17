use bio_types::strand::Strand;
use derive_more::{Add, AddAssign};

use crate::core::summary::{MismatchesSummary, ROISummary};

use super::ROIBasedStat;

#[derive(Copy, Clone, Default, Add, AddAssign)]
pub struct EditingIndex {
    accumulator: MismatchesSummary,
}

impl ROIBasedStat for EditingIndex {
    fn process(&mut self, summary: &ROISummary) {
        match summary.strand {
            Strand::Forward => {
                self.accumulator += summary.mismatches;
            }
            Strand::Reverse => {
                self.accumulator += summary.mismatches.complementary();
            }
            Strand::Unknown => {}
        }
    }

    fn combine(stats: &[Self]) -> Self {
        stats.into_iter().fold(Self::default(), |a, b| a + *b)
    }

    fn header() -> &'static str {
        "A->A\t\
        T->T\t\
        G->G\t\
        C->C\t\
        A->T\t\
        T->A\t\
        A->G\t\
        T->C\t\
        A->C\t\
        T->G\t\
        G->C\t\
        C->G\t\
        G->A\t\
        C->T\t\
        G->T\t\
        C->A"
    }

    fn row(&self) -> String {
        let res = &self.accumulator;
        format!(
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            res.A.A as f32 / res.A.coverage() as f32,
            res.T.T as f32 / res.T.coverage() as f32,
            res.G.G as f32 / res.G.coverage() as f32,
            res.C.C as f32 / res.C.coverage() as f32,
            res.A.T as f32 / res.A.coverage() as f32,
            res.T.A as f32 / res.T.coverage() as f32,
            res.A.G as f32 / res.A.coverage() as f32,
            res.T.C as f32 / res.T.coverage() as f32,
            res.A.C as f32 / res.A.coverage() as f32,
            res.T.G as f32 / res.T.coverage() as f32,
            res.G.C as f32 / res.G.coverage() as f32,
            res.C.G as f32 / res.C.coverage() as f32,
            res.G.A as f32 / res.G.coverage() as f32,
            res.C.T as f32 / res.C.coverage() as f32,
            res.G.T as f32 / res.G.coverage() as f32,
            res.C.A as f32 / res.C.coverage() as f32
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::core::counting::NucCounts;
    use bio_types::genome::Interval;

    #[test]
    fn header() {
        assert_eq!(
            EditingIndex::header(),
            "A->A\tT->T\tG->G\tC->C\tA->T\tT->A\tA->G\tT->C\tA->C\tT->G\tG->C\tC->G\tG->A\tC->T\tG->T\tC->A"
        )
    }

    #[test]
    fn row() {
        let mut dummy = EditingIndex::default();

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

        dummy.process(&ROISummary {
            interval: Interval::new("".into(), 1..2),
            strand: Strand::Forward,
            name: "".into(),
            sequenced: NucCounts::new(1, 1, 1, 1),
            mismatches: forward,
        });

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

        dummy.process(&ROISummary {
            interval: Interval::new("".into(), 1..2),
            strand: Strand::Reverse,
            name: "".into(),
            sequenced: NucCounts::new(1, 1, 1, 1),
            mismatches: reverse,
        });

        assert_eq!(
            dummy.row(),
            format!(
                "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                (10 + 17) as f32 / (19 + 62) as f32,
                (16 + 2) as f32 / (58 + 14) as f32,
                (11 + 7) as f32 / (42 + 30) as f32,
                (6 + 12) as f32 / (36 + 46) as f32,
                (4 + 14) as f32 / (19 + 62) as f32,
                (13 + 5) as f32 / (58 + 14) as f32,
                (3 + 15) as f32 / (19 + 62) as f32,
                (14 + 4) as f32 / (58 + 14) as f32,
                (2 + 16) as f32 / (19 + 62) as f32,
                (15 + 3) as f32 / (58 + 14) as f32,
                (10 + 8) as f32 / (42 + 30) as f32,
                (7 + 11) as f32 / (36 + 46) as f32,
                (9 + 9) as f32 / (42 + 30) as f32,
                (8 + 10) as f32 / (36 + 46) as f32,
                (12 + 6) as f32 / (42 + 30) as f32,
                (15 + 13) as f32 / (36 + 46) as f32,
            )
        )
    }
}
