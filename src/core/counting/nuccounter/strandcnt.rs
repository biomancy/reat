use std::cell::RefCell;
use std::cmp::{max, min};
use std::marker::PhantomData;
use std::ops::{Index, Range};

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{ReqStrand, Strand};
use itertools::{chain, zip, Itertools};
use rust_htslib::bam::record::Cigar;

use crate::core::counting::buffers::IntervalCounts;
use crate::core::counting::nuccounter::{CountingResult, NucCounter};
use crate::core::counting::NucCounts;
use crate::core::filtering::reads::ReadsFilter;
use crate::core::read::AlignedRead;
use crate::core::stranding::deduct::StrandDeductor;
use crate::core::workload::ROIWorkload;

use super::super::buffers::CountsBuffer;

#[derive(Clone)]
pub struct StrandedNucCounter<R: AlignedRead, Deductor: StrandDeductor<R>, Counter: NucCounter<R>> {
    forward: Counter,
    reverse: Counter,
    sequenced: Vec<NucCounts>,
    deductor: Deductor,
    phantom: PhantomData<fn() -> R>,
}

impl<R: AlignedRead, Deductor: StrandDeductor<R>, Counter: NucCounter<R>> StrandedNucCounter<R, Deductor, Counter> {
    pub fn new(forward: Counter, reverse: Counter, deductor: Deductor) -> Self {
        debug_assert_eq!(forward.rois(), reverse.rois());
        debug_assert_eq!(forward.interval(), reverse.interval());
        Self { forward, reverse, sequenced: vec![], deductor, phantom: Default::default() }
    }
}

impl<R: AlignedRead, Deductor: StrandDeductor<R>, Counter: NucCounter<R>> NucCounter<R>
    for StrandedNucCounter<R, Deductor, Counter>
{
    #[inline]
    fn interval(&self) -> &Interval {
        debug_assert!(self.forward.interval() == self.reverse.interval());
        self.forward.interval()
    }

    #[inline]
    fn rois(&self) -> &[Interval] {
        debug_assert_eq!(self.forward.rois(), self.reverse.rois());
        &self.forward.rois()
    }

    fn counted(&mut self) -> &[NucCounts] {
        let (forward, reverse) = (self.forward.counted(), self.reverse.counted());
        debug_assert_eq!(forward.len(), reverse.len());
        self.sequenced = zip(forward, reverse).map(|(a, b)| *a + *b).collect();
        &self.sequenced
    }

    #[inline]
    fn mapped(&self) -> u32 {
        self.forward.mapped() + self.reverse.mapped()
    }

    fn count(&mut self, read: &mut R) {
        match self.deductor.deduce(&read) {
            ReqStrand::Forward => self.forward.count(read),
            ReqStrand::Reverse => self.reverse.count(read),
        }
    }

    fn results(&self) -> Vec<CountingResult> {
        let (forward, reverse) = (self.forward.results(), self.reverse.results());
        let mut result = Vec::with_capacity(forward.len() + reverse.len());
        for (mut f, mut r) in zip(forward, reverse) {
            debug_assert_eq!(f.roi, r.roi);
            debug_assert_eq!(f.name, r.name);
            f.strand = Strand::Forward;
            result.push(f);
            r.strand = Strand::Reverse;
            result.push(r);
        }
        result
    }

    #[inline]
    fn reset(&mut self, info: ROIWorkload) {
        self.forward.reset(info.clone());
        self.reverse.reset(info);
        self.sequenced.clear();
    }
}

#[cfg(test)]
mod tests {
    use mockall::{predicate::eq, Sequence};

    use crate::core::counting::buffers::{MockCountsBuffer, RawCounts};
    use crate::core::counting::nuccounter::MockNucCounter;
    use crate::core::read::MockRead;
    use crate::core::stranding::deduct::MockStrandDeductor;

    use super::*;

    fn counter() -> MockNucCounter<MockRead> {
        let mut counter = MockNucCounter::<MockRead>::new();
        counter.expect_rois().return_const(vec![Interval::new("".into(), 0..1), Interval::new("".into(), 0..2)]);
        counter.expect_interval().return_const(Interval::new("".into(), 0..5));
        counter
    }

    #[test]
    fn mapped() {
        let mut forward = counter();
        forward.expect_mapped().return_const(3u32);

        let mut reverse = counter();
        reverse.expect_mapped().return_const(123u32);

        let dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
        debug_assert_eq!(dummy.mapped(), 126);
    }

    #[test]
    fn delegation() {
        let dummy = StrandedNucCounter::new(counter(), counter(), MockStrandDeductor::new());

        let reference = counter();
        debug_assert_eq!(dummy.rois(), reference.rois());
        debug_assert_eq!(dummy.interval(), reference.interval());
    }

    #[test]
    fn count() {
        // Forward
        let mut seq = Sequence::new();

        let mut deductor = MockStrandDeductor::new();
        let mut forward = counter();
        let mut reverse = counter();

        deductor.expect_deduce().once().return_const(ReqStrand::Forward).in_sequence(&mut seq);
        forward.expect_count().once().return_const(()).in_sequence(&mut seq);

        deductor.expect_deduce().once().return_const(ReqStrand::Reverse).in_sequence(&mut seq);
        reverse.expect_count().once().return_const(()).in_sequence(&mut seq);

        let mut dummy = StrandedNucCounter::new(forward, reverse, deductor);

        let mut read = MockRead::new();
        dummy.count(&mut read);
        dummy.count(&mut read);
    }

    #[test]
    fn reset() {
        let mut forward = counter();
        forward.expect_reset().once().return_const(());

        let mut reverse = counter();
        reverse.expect_reset().once().return_const(());

        let mut dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
        dummy.sequenced.extend(&[NucCounts::A(123), NucCounts::C(1), NucCounts::G(1233)]);

        dummy.reset(ROIWorkload::default());
        assert!(dummy.sequenced.is_empty());
    }

    #[test]
    fn counted() {
        let mut forward = counter();
        forward.expect_counted().once().return_const(vec![NucCounts::A(12), NucCounts::G(3), NucCounts::T(2)]);

        let mut reverse = counter();
        reverse.expect_counted().once().return_const(vec![NucCounts::A(3), NucCounts::G(7), NucCounts::G(1)]);

        let mut dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
        let counted = dummy.counted();

        assert_eq!(counted, &[NucCounts::A(15), NucCounts::G(10), NucCounts { A: 0, C: 0, G: 1, T: 2 }]);
    }

    const NUCLEOTIDES_1_1: &[NucCounts] = &[NucCounts::A(1), NucCounts::G(3)];
    const NUCLEOTIDES_1_2: &[NucCounts] = &[NucCounts::G(33), NucCounts::C(31)];
    const NUCLEOTIDES_2_1: &[NucCounts] = &[NucCounts::A(1), NucCounts::T(3), NucCounts::T(3)];
    const NUCLEOTIDES_2_2: &[NucCounts] = &[NucCounts::C(1), NucCounts::T(3), NucCounts::G(3)];

    #[test]
    fn results() {
        // Leak it intentionally
        let interval1 = Box::leak(Box::new(Interval::new("3".into(), 0..124)));
        let interval2 = Box::leak(Box::new(Interval::new("4".into(), 10..124)));

        let mut fwd = vec![
            CountingResult {
                name: "1",
                strand: Strand::Unknown,
                roi: interval1,
                cnts: RawCounts { nuc: NUCLEOTIDES_1_1, coverage: 5 },
            },
            CountingResult {
                name: "2",
                strand: Strand::Unknown,
                roi: interval2,
                cnts: RawCounts { nuc: NUCLEOTIDES_2_2, coverage: 5123 },
            },
        ];
        let mut forward = counter();
        forward.expect_results().once().return_const(fwd.clone());

        let mut rev = vec![
            CountingResult {
                name: "1",
                strand: Strand::Forward,
                roi: interval1,
                cnts: RawCounts { nuc: NUCLEOTIDES_1_2, coverage: 55 },
            },
            CountingResult {
                name: "2",
                strand: Strand::Reverse,
                roi: interval2,
                cnts: RawCounts { nuc: NUCLEOTIDES_2_1, coverage: 15 },
            },
        ];
        let mut reverse = counter();
        reverse.expect_results().once().return_const(rev.clone());

        let mut dummy = StrandedNucCounter::new(forward, reverse, MockStrandDeductor::new());
        let mut counted = dummy.results();
        counted.sort_by_key(|x| x.name);

        fwd.extend(rev);
        let mut expected = fwd
            .into_iter()
            .zip(&[Strand::Forward, Strand::Forward, Strand::Reverse, Strand::Reverse])
            .map(|(mut x, strand)| {
                x.strand = *strand;
                x
            })
            .collect_vec();
        expected.sort_by_key(|x| x.name);
        assert_eq!(counted, expected);
    }
}
