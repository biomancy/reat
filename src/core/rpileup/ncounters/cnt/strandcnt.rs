use std::cell::RefCell;
use std::marker::PhantomData;
use std::ops::Range;
use std::rc::Rc;
use std::sync::Mutex;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::{zip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::BatchedMismatches;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, NucCounter};
use crate::core::rpileup::ReadsCollider;
use crate::core::stranding::deduct::StrandDeductor;

pub struct StrandedNucCounter<'a, R, Deductor, InnerNucCounter>
where
    R: AlignedRead,
    Deductor: StrandDeductor<R>,
    InnerNucCounter: NucCounter<'a, R>,
    InnerNucCounter::ColliderResult: AggregatedNucCounts<'a>,
{
    forward: InnerNucCounter,
    reverse: InnerNucCounter,
    deductor: Deductor,
    phantom: PhantomData<(&'a R)>,
}

impl<'a, R, Deductor, InnerNucCounter> ReadsCollider<'a, R> for StrandedNucCounter<'a, R, Deductor, InnerNucCounter>
where
    R: AlignedRead,
    Deductor: StrandDeductor<R>,
    InnerNucCounter: NucCounter<'a, R>,
    InnerNucCounter::Workload: Clone,
    InnerNucCounter::ColliderResult: AggregatedNucCounts<'a>,
    <<InnerNucCounter as ReadsCollider<'a, R>>::ColliderResult as AggregatedNucCounts<'a>>::ItemInfo: PartialEq,
{
    type ColliderResult = InnerNucCounter::ColliderResult;
    type Workload = InnerNucCounter::Workload;

    fn reset(&mut self, info: Self::Workload) {
        self.forward.reset(info.clone());
        self.reverse.reset(info);
    }

    #[inline]
    fn collide(&mut self, read: &R) {
        match self.deductor.deduce(read) {
            ReqStrand::Forward => self.forward.collide(read),
            ReqStrand::Reverse => self.reverse.collide(read),
        }
    }

    fn finalize(&mut self) {
        self.forward.finalize();
        self.reverse.finalize();
    }

    fn result(&'a self) -> Self::ColliderResult {
        let (mut fwd, mut rev) = (self.forward.result(), self.reverse.result());

        for (f, r) in zip(fwd.items_mut(), rev.items_mut()) {
            debug_assert!(f.unstranded.is_some() && r.unstranded.is_none());
            debug_assert!(f.forward.is_none() && f.reverse.is_none() && r.forward.is_none() && r.reverse.is_none());
            debug_assert!(f.info == r.info);

            f.forward = f.unstranded;
            f.reverse = r.unstranded;
            f.unstranded = None;
        }
        fwd
    }
}

// impl<'a, R, Deductor, InnerNucCounter> NucCounter<'a, R> for StrandedNucCounter<'a, R, Deductor, InnerNucCounter>
// where
//     R: AlignedRead,
//     Deductor: StrandDeductor<R>,
//     InnerNucCounter: NucCounter<'a, R>,
//     InnerNucCounter::Workload: Clone,
//     InnerNucCounter::ColliderResult: AggregatedNucCounts<'a>,
//     <<InnerNucCounter as ReadsCollider<'a, R>>::ColliderResult as AggregatedNucCounts<'a>>::ItemInfo: PartialEq,
// {
// }

// Fetch results
// let ((finter, forward), (frev, reverse)) = (self.forward.result().dissolve(), self.reverse.result().dissolve());
// debug_assert_eq!(finter, frev);
//
// // Calculate shared counts once to avoid repeated memory allocations
// let cnts = Rc::new(
//     zip(&forward, &reverse)
//         .map(|(f, r)| zip(f.ncounts(), r.ncounts()).map(|(c1, c2)| *c1 + *c2))
//         .flatten()
//         .collect_vec(),
// );
//
// let mut results = Vec::with_capacity(forward.len());
// let mut start = 0;
// for (f, r) in zip(forward, reverse) {
//     debug_assert_eq!(f.contig(), r.contig());
//     debug_assert_eq!(f.range(), r.range());
//
//     let end = start + f.ncounts().len();
//     results.push(StrandedNucCounts {
//         forward: f,
//         reverse: r,
//         idx: start..end,
//         ncounts: cnts.clone(),
//         marker: Default::default(),
//     });
//     start = end;
// }

// #[cfg(test)]
// mod tests {
//     use mockall::Sequence;
//
//     use crate::core::counting::buffers::RawCounts;
//     use crate::core::counting::nuccounter::MockNucCounter;
//     use crate::core::read::MockRead;
//     use crate::core::stranding::deduct::MockStrandDeductor;
//
//     use super::*;
//     use itertools::Itertools;
//
//     fn counter() -> MockNucCounter<MockRead> {
//         let mut counter = MockNucCounter::<MockRead>::new();
//         counter.expect_rois().return_const(vec![Interval::new("".into(), 0..1), Interval::new("".into(), 0..2)]);
//         counter.expect_interval().return_const(Interval::new("".into(), 0..5));
//         counter
//     }
//
//     #[test]
//     fn mapped() {
//         let mut forward = counter();
//         forward.expect_mapped().return_const(3u32);
//
//         let mut reverse = counter();
//         reverse.expect_mapped().return_const(123u32);
//
//         let dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
//         debug_assert_eq!(dummy.mapped(), 126);
//     }
//
//     #[test]
//     fn delegation() {
//         let dummy = StrandedNucCounter::new(counter(), counter(), MockStrandDeductor::new());
//
//         let reference = counter();
//         debug_assert_eq!(dummy.rois(), reference.rois());
//         debug_assert_eq!(dummy.interval(), reference.interval());
//     }
//
//     #[test]
//     fn count() {
//         // Forward
//         let mut seq = Sequence::new();
//
//         let mut deductor = MockStrandDeductor::new();
//         let mut forward = counter();
//         let mut reverse = counter();
//
//         deductor.expect_deduce().once().return_const(ReqStrand::Forward).in_sequence(&mut seq);
//         forward.expect_count().once().return_const(()).in_sequence(&mut seq);
//
//         deductor.expect_deduce().once().return_const(ReqStrand::Reverse).in_sequence(&mut seq);
//         reverse.expect_count().once().return_const(()).in_sequence(&mut seq);
//
//         let mut dummy = StrandedNucCounter::new(forward, reverse, deductor);
//
//         let mut read = MockRead::new();
//         dummy.count(&mut read);
//         dummy.count(&mut read);
//     }
//
//     #[test]
//     fn reset() {
//         let mut forward = counter();
//         forward.expect_reset().once().return_const(());
//
//         let mut reverse = counter();
//         reverse.expect_reset().once().return_const(());
//
//         let mut dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
//         dummy.sequenced.extend(&[NucCounts::A(123), NucCounts::C(1), NucCounts::G(1233)]);
//
//         dummy.reset(ROIWorkload::default());
//         assert!(dummy.sequenced.is_empty());
//     }
//
//     #[test]
//     fn counted() {
//         let mut forward = counter();
//         forward.expect_counted().once().return_const(vec![NucCounts::A(12), NucCounts::G(3), NucCounts::T(2)]);
//
//         let mut reverse = counter();
//         reverse.expect_counted().once().return_const(vec![NucCounts::A(3), NucCounts::G(7), NucCounts::G(1)]);
//
//         let mut dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
//         let counted = dummy.counted();
//
//         assert_eq!(counted, &[NucCounts::A(15), NucCounts::G(10), NucCounts { A: 0, C: 0, G: 1, T: 2 }]);
//     }
//
//     const NUCLEOTIDES_1_1: &[NucCounts] = &[NucCounts::A(1), NucCounts::G(3)];
//     const NUCLEOTIDES_1_2: &[NucCounts] = &[NucCounts::G(33), NucCounts::C(31)];
//     const NUCLEOTIDES_2_1: &[NucCounts] = &[NucCounts::A(1), NucCounts::T(3), NucCounts::T(3)];
//     const NUCLEOTIDES_2_2: &[NucCounts] = &[NucCounts::C(1), NucCounts::T(3), NucCounts::G(3)];
//
//     #[test]
//     fn results() {
//         // Leak it intentionally
//         let interval1 = Box::leak(Box::new(Interval::new("3".into(), 0..124)));
//         let interval2 = Box::leak(Box::new(Interval::new("4".into(), 10..124)));
//
//         let mut fwd = vec![
//             CountingResult {
//                 name: "1",
//                 strand: Strand::Unknown,
//                 roi: interval1,
//                 cnts: RawCounts { nuc: NUCLEOTIDES_1_1, coverage: 5 },
//             },
//             CountingResult {
//                 name: "2",
//                 strand: Strand::Unknown,
//                 roi: interval2,
//                 cnts: RawCounts { nuc: NUCLEOTIDES_2_2, coverage: 5123 },
//             },
//         ];
//         let mut forward = counter();
//         forward.expect_results().once().return_const(fwd.clone());
//
//         let rev = vec![
//             CountingResult {
//                 name: "1",
//                 strand: Strand::Forward,
//                 roi: interval1,
//                 cnts: RawCounts { nuc: NUCLEOTIDES_1_2, coverage: 55 },
//             },
//             CountingResult {
//                 name: "2",
//                 strand: Strand::Reverse,
//                 roi: interval2,
//                 cnts: RawCounts { nuc: NUCLEOTIDES_2_1, coverage: 15 },
//             },
//         ];
//         let mut reverse = counter();
//         reverse.expect_results().once().return_const(rev.clone());
//
//         let dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
//         let mut counted = dummy.results();
//         counted.sort_by_key(|x| x.name);
//
//         fwd.extend(rev);
//         let mut expected = fwd
//             .into_iter()
//             .zip(&[Strand::Forward, Strand::Forward, Strand::Reverse, Strand::Reverse])
//             .map(|(mut x, strand)| {
//                 x.strand = *strand;
//                 x
//             })
//             .collect_vec();
//         expected.sort_by_key(|x| x.name);
//         assert_eq!(counted, expected);
//     }
// }
