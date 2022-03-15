use std::any::Any;

use bio_types::strand::Strand;
use derive_more::{Add, AddAssign};

use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::core::hooks::stats::EditingStat;
use crate::core::hooks::stats::EditingStatType;
use crate::core::hooks::Hook;
use crate::core::mismatches::roi::{NucMismatches, ROIMismatchesVec};
use crate::core::mismatches::{Batch};

#[derive(Copy, Clone, Default, Add, AddAssign)]
pub struct ROIEditingIndex {
    accumulator: NucMismatches,
    unstranded: usize,
}

impl ROIEditingIndex {
    fn process(&mut self, x: &ROIMismatchesVec, strand: Strand) {
        let iter = x.data.mismatches.iter();

        match strand {
            Strand::Forward => {
                for x in iter {
                    self.accumulator += *x;
                }
            }
            Strand::Reverse => {
                for x in iter {
                    self.accumulator += x.complementary();
                }
            }
            Strand::Unknown => {
                self.unstranded += 1;
            }
        }
    }

    pub fn collapse(items: Vec<Box<dyn Any>>) -> Self {
        items.into_iter().map(|x| x.downcast::<Self>().unwrap()).fold(Self::default(), |a, b| a + *b)
    }
}

impl Serialize for ROIEditingIndex {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let res = &self.accumulator;

        let mut state = serializer.serialize_struct("ROIEditingIndex", 16)?;
        state.serialize_field("A->A", &(res.A.A as f32 / res.A.coverage() as f32))?;
        state.serialize_field("T->T", &(res.T.T as f32 / res.T.coverage() as f32))?;
        state.serialize_field("G->G", &(res.G.G as f32 / res.G.coverage() as f32))?;
        state.serialize_field("C->C", &(res.C.C as f32 / res.C.coverage() as f32))?;
        state.serialize_field("A->T", &(res.A.T as f32 / res.A.coverage() as f32))?;
        state.serialize_field("T->A", &(res.T.A as f32 / res.T.coverage() as f32))?;
        state.serialize_field("A->G", &(res.A.G as f32 / res.A.coverage() as f32))?;
        state.serialize_field("T->C", &(res.T.C as f32 / res.T.coverage() as f32))?;
        state.serialize_field("A->C", &(res.A.C as f32 / res.A.coverage() as f32))?;
        state.serialize_field("T->G", &(res.T.G as f32 / res.T.coverage() as f32))?;
        state.serialize_field("G->C", &(res.G.C as f32 / res.G.coverage() as f32))?;
        state.serialize_field("C->G", &(res.C.G as f32 / res.C.coverage() as f32))?;
        state.serialize_field("G->A", &(res.G.A as f32 / res.G.coverage() as f32))?;
        state.serialize_field("C->T", &(res.C.T as f32 / res.C.coverage() as f32))?;
        state.serialize_field("G->T", &(res.G.T as f32 / res.G.coverage() as f32))?;
        state.serialize_field("C->A", &(res.C.A as f32 / res.C.coverage() as f32))?;
        state.end()
    }
}

impl Hook<ROIMismatchesVec> for ROIEditingIndex {
    fn on_finish(&mut self, mismatches: &mut Batch<ROIMismatchesVec>) {
        for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
            self.process(&mismatches.retained[strand], strand);
            self.process(&mismatches.items[strand], strand);
        }
    }
}

impl EditingStat<ROIMismatchesVec> for ROIEditingIndex {
    fn into_any(self: Box<Self>) -> (EditingStatType, Box<dyn Any>) {
        (EditingStatType::ROIEditingIndex, self)
    }
}

// #[cfg(test)]
// mod test {
//     use crate::core::mismatches::roi::MockBatchedROIMismatches;
//
//     use super::*;
//
//     #[test]
//     fn row() {
// let factory = |roi, strand, coverage, sequence, mismatches| {
//     let mut dummy = MockROIMismatches::new();
//     dummy.expect_roi().return_const(roi);
//     dummy.expect_strand().return_const(strand);
//     dummy.expect_coverage().return_const(coverage);
//     dummy.expect_sequence().return_const(sequence);
//     dummy.expect_mismatches().return_const(mismatches);
//     dummy
// };
//
// let mut dummy = ROIEditingIndex::default();
//
// let mut forward = MismatchesSummary::zeros();
// // coverage - 19
// forward.A.A = 10;
// forward.A.C = 2;
// forward.A.G = 3;
// forward.A.T = 4;
//
// // coverage - 36
// forward.C.A = 15;
// forward.C.C = 6;
// forward.C.G = 7;
// forward.C.T = 8;
//
// // coverage - 42
// forward.G.A = 9;
// forward.G.C = 10;
// forward.G.G = 11;
// forward.G.T = 12;
//
// // coverage - 58
// forward.T.A = 13;
// forward.T.C = 14;
// forward.T.G = 15;
// forward.T.T = 16;
//
// let forward = factory(
//     ROI::new(Interval::new("".into(), 1..2), "".into(), vec![1..2]),
//     Strand::Forward,
//     3u32,
//     NucCounts::new(1, 1, 1, 1),
//     forward,
// );
//
// let mut reverse = MismatchesSummary::zeros();
// // coverage - 14
// reverse.A.A = 2;
// reverse.A.C = 3;
// reverse.A.G = 4;
// reverse.A.T = 5;
//
// // coverage - 30
// reverse.C.A = 6;
// reverse.C.C = 7;
// reverse.C.G = 8;
// reverse.C.T = 9;
//
// // coverage - 46
// reverse.G.A = 10;
// reverse.G.C = 11;
// reverse.G.G = 12;
// reverse.G.T = 13;
//
// // coverage - 62
// reverse.T.A = 14;
// reverse.T.C = 15;
// reverse.T.G = 16;
// reverse.T.T = 17;
//
// let reverse = factory(
//     ROI::new(Interval::new("".into(), 1..2), "".into(), vec![1..2]),
//     Strand::Reverse,
//     9u32,
//     NucCounts::new(1, 1, 1, 1),
//     reverse,
// );
//
// let prehook = vec![forward, reverse];
// let posthook = dummy.hook(prehook);
//
// assert_eq!(
//     dummy.row(),
//     vec!(
//         format!("{:.6}", (10 + 17) as f32 / (19 + 62) as f32),
//         format!("{:.6}", (16 + 2) as f32 / (58 + 14) as f32),
//         format!("{:.6}", (11 + 7) as f32 / (42 + 30) as f32),
//         format!("{:.6}", (6 + 12) as f32 / (36 + 46) as f32),
//         format!("{:.6}", (4 + 14) as f32 / (19 + 62) as f32),
//         format!("{:.6}", (13 + 5) as f32 / (58 + 14) as f32),
//         format!("{:.6}", (3 + 15) as f32 / (19 + 62) as f32),
//         format!("{:.6}", (14 + 4) as f32 / (58 + 14) as f32),
//         format!("{:.6}", (2 + 16) as f32 / (19 + 62) as f32),
//         format!("{:.6}", (15 + 3) as f32 / (58 + 14) as f32),
//         format!("{:.6}", (10 + 8) as f32 / (42 + 30) as f32),
//         format!("{:.6}", (7 + 11) as f32 / (36 + 46) as f32),
//         format!("{:.6}", (9 + 9) as f32 / (42 + 30) as f32),
//         format!("{:.6}", (8 + 10) as f32 / (36 + 46) as f32),
//         format!("{:.6}", (12 + 6) as f32 / (42 + 30) as f32),
//         format!("{:.6}", (15 + 13) as f32 / (36 + 46) as f32),
//     )
// )
// }
// }
