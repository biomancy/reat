use std::cmp::min;
use std::marker::PhantomData;

use bio_types::genome::{AbstractInterval, Interval};
use rust_htslib::bam::record::Cigar;

use crate::core::counting::LocusCounts;
use crate::core::filtering::reads::ReadsFilter;
use crate::core::read::AlignedRead;

use super::super::buffers::CountsBuffer;
use super::{NucCounter, NucCounterContent};

#[derive(Clone)]
pub struct BaseNucCounter<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>> {
    filter: Filter,
    buffer: Buffer,
    roi: Interval,
    reads_counted: usize,
    phantom: PhantomData<R>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>> BaseNucCounter<R, Filter, Buffer> {
    #[inline]
    fn is_record_ok(&self, record: &mut R) -> bool {
        self.filter.is_read_ok(record) && record.contig() == self.roi.contig()
    }

    pub fn new(filter: Filter, mut buffer: Buffer, roi: Interval) -> Self {
        buffer.reset((roi.range().end - roi.range().start) as u32);
        BaseNucCounter { filter, buffer, roi, reads_counted: 0, phantom: Default::default() }
    }

    fn implprocess(&mut self, read: &mut R) -> &[LocusCounts] {
        let counts = self.buffer.buffer_for(read);
        let sequence = read.seq();

        let (mut roipos, mut seqpos) = (read.pos() - self.roi.range().start as i64, 0usize);
        let roisize = (self.roi.range().end - self.roi.range().start) as i64;

        for cigar in read.cigar().iter() {
            if roipos >= roisize {
                break;
            }
            match cigar {
                Cigar::Match(ops) | Cigar::Equal(ops) | Cigar::Diff(ops) => {
                    // fast-forward when possible
                    let skip = if roipos < 0 {
                        let skip = min(roipos.abs() as u32, *ops);
                        roipos += skip as i64;
                        seqpos += skip as usize;
                        skip
                    } else {
                        0
                    };
                    for _ in skip..*ops {
                        if roipos >= roisize {
                            break;
                        } else if self.filter.is_base_ok(read, seqpos) {
                            debug_assert!(roipos >= 0);
                            // From the SAM specification: No assumptions can be made on the letter cases
                            match sequence[seqpos as usize] {
                                b'A' | b'a' => counts[roipos as usize].A += 1,
                                b'T' | b't' => counts[roipos as usize].T += 1,
                                b'G' | b'g' => counts[roipos as usize].G += 1,
                                b'C' | b'c' => counts[roipos as usize].C += 1,
                                _ => {}
                            }
                        }
                        roipos += 1;
                        seqpos += 1;
                    }
                }
                Cigar::Del(ops) | Cigar::RefSkip(ops) => {
                    roipos += *ops as i64;
                }
                Cigar::SoftClip(ops) | Cigar::Ins(ops) => {
                    seqpos += *ops as usize;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        counts
    }
}

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>> NucCounter<R>
    for BaseNucCounter<R, Filter, Buffer>
{
    #[inline]
    fn roi(&self) -> &Interval {
        &self.roi
    }

    #[inline]
    fn process(&mut self, read: &mut R) {
        if !self.is_record_ok(read) {
            return;
        }
        self.implprocess(read);
        self.reads_counted += 1;
    }

    fn reads_counted(&self) -> usize {
        self.reads_counted
    }

    #[inline]
    fn reset(&mut self, roi: Interval) {
        let size = (roi.range().end - roi.range().start) as u32;
        self.roi = roi;
        self.buffer.reset(size);
        self.reads_counted = 0;
    }

    #[inline]
    fn content(&self) -> NucCounterContent {
        NucCounterContent { interval: self.roi.clone(), counts: self.buffer.content() }
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use mockall::predicate;
    use rust_htslib::bam::record::Cigar::*;
    use rust_htslib::bam::record::CigarString;

    use crate::core::counting::{CountsBufferContent, LocusCounts};
    use crate::core::filtering::reads::MockReadsFilter;
    use crate::core::read::MockRead;

    use super::super::super::buffers::MockCountsBuffer;
    use super::*;

    fn _implprocess(
        roi: Range<u64>,
        pos: i64,
        seq: &str,
        okbase: Vec<bool>,
        cigar: Vec<Cigar>,
        expected: &[LocusCounts],
    ) {
        // vec![Cigar::Match(100), Cigar::SoftClip(10)]
        let roi = Interval::new("".into(), roi);
        let roisize = roi.range().end - roi.range().start;

        let mut buffer = MockCountsBuffer::<MockRead>::new();
        buffer.expect_reset().once().return_const(());
        buffer.expect_buffer_for().once().returning(move |_| vec![LocusCounts::zeros()].repeat(roisize as usize));

        let mut filter = MockReadsFilter::new();
        for isok in okbase {
            filter.expect_is_base_ok().once().return_const(isok);
        }

        let mut counter = BaseNucCounter::new(filter, buffer, roi);

        let mut read = MockRead::new();
        read.expect_pos().return_const(pos);
        read.expect_cigar().return_once(move || CigarString(cigar).into_view(pos));
        let seq = String::from(seq);
        read.expect_seq().returning(move || seq.as_bytes().to_vec());

        let result = counter.implprocess(&mut read);
        assert_eq!(result.len(), expected.len());
        assert!(result.iter().zip(expected).all(|(x, y)| x == y));
    }

    #[test]
    #[allow(non_snake_case)]
    fn implprocess() {
        let A = || LocusCounts::new(1, 0, 0, 0);
        let C = || LocusCounts::new(0, 1, 0, 0);
        let G = || LocusCounts::new(0, 0, 1, 0);
        let T = || LocusCounts::new(0, 0, 0, 1);
        let Z = || LocusCounts::zeros();

        // Query consuming operations only
        let M = |x| Match(x);
        let X = |x| Diff(x);
        let E = |x| Equal(x);
        for op in [M, X, E] {
            // complete overlap with the region
            _implprocess(0..4, 0, "ACGT", vec![true, true, true, true], vec![op(4)], &vec![A(), C(), G(), T()]);
            // inside region
            _implprocess(0..4, 1, "AC", vec![true, true], vec![op(2)], &vec![Z(), A(), C(), Z()]);
            // completely out of the region
            _implprocess(0..4, 4, "ACGT", vec![], vec![op(4)], &vec![Z(), Z(), Z(), Z()]);
            // end out of the region
            _implprocess(0..4, 2, "ACGT", vec![true, true], vec![op(4)], &vec![Z(), Z(), A(), C()]);
            // start out of the region
            _implprocess(0..4, -1, "ACGT", vec![true, true, true], vec![op(4)], &vec![C(), G(), T(), Z()]);
        }

        // No-ops + reference/query consuming operations only
        let D = |x| Del(x);
        let N = |x| RefSkip(x);
        let H = |x| HardClip(x);
        let P = |x| Pad(x);
        let S = |x| SoftClip(x);
        let I = |x| Ins(x);

        let empty = vec![Z()].repeat(4);
        for op in [D, N, H, P, S, I] {
            // complete overlap with the region
            _implprocess(0..4, 0, "ACGT", vec![], vec![op(4)], &empty);
            // inside region
            _implprocess(0..4, 1, "AC", vec![], vec![op(2)], &empty);
            // completely out of the region
            _implprocess(0..4, 4, "ACGT", vec![], vec![op(4)], &empty);
            // end out of the region
            _implprocess(0..4, 2, "ACGT", vec![], vec![op(4)], &empty);
            // start out of the region
            _implprocess(0..4, -1, "ACGT", vec![], vec![op(4)], &empty);
        }

        // Complex queries
        _implprocess(2..5, 0, "AGC", vec![true, true], vec![D(2), M(1), N(1), M(1), N(1), M(1)], &vec![A(), Z(), G()]);
        _implprocess(2..5, 0, "AGC", vec![true], vec![D(4), M(1), N(3), M(2)], &vec![Z(), Z(), A()]);
        _implprocess(2..5, 0, "AGC", vec![false], vec![M(3)], &vec![Z(), Z(), Z()]);
        _implprocess(2..5, 0, "AGC", vec![true, true], vec![M(1), N(2), M(2)], &vec![Z(), G(), C()]);
        _implprocess(1..5, 0, "NNN", vec![true, true], vec![M(1), N(2), M(2)], &vec![Z(), Z(), Z(), Z()]);
    }

    #[test]
    fn reset() {
        let mut buffer = MockCountsBuffer::<MockRead>::new();
        buffer.expect_reset().with(predicate::eq(100u32)).once().return_const(());
        buffer.expect_reset().with(predicate::eq(1000u32)).once().return_const(());
        let mut dummy = BaseNucCounter::new(MockReadsFilter::new(), buffer, Interval::new("".into(), 100..200));

        let newroi = Interval::new("".to_string(), 1000..2000);
        dummy.reset(newroi);
    }

    #[test]
    fn content() {
        let mut buffer = MockCountsBuffer::<MockRead>::new();
        buffer.expect_reset().once().return_const(());
        buffer.expect_content().once().return_const(CountsBufferContent {
            forward: None,
            reverse: None,
            unstranded: None,
        });
        let dummy = BaseNucCounter::new(MockReadsFilter::new(), buffer, Interval::new("".into(), 1..2));

        let content = dummy.content();
        assert_eq!(content.interval, dummy.roi);
        assert_eq!(content.counts.forward, None);
        assert_eq!(content.counts.reverse, None);
        assert_eq!(content.counts.unstranded, None);
    }

    #[test]
    fn is_record_ok() {
        let contig = "".to_string();
        let wrong_contig = "!".to_string();
        let roi = Interval::new(contig.clone(), 100..200);
        for (isok, ctg, result) in [
            (true, &contig, true),
            (false, &contig, false),
            (true, &wrong_contig, false),
            (false, &wrong_contig, false),
        ] {
            let mut filter = MockReadsFilter::new();
            filter.expect_is_read_ok().once().return_const(isok);

            let mut buffer = MockCountsBuffer::<MockRead>::new();
            buffer.expect_reset().once().return_const(());
            let dummy = BaseNucCounter::new(filter, buffer, roi.clone());

            let mut read = MockRead::new();
            read.expect_contig().return_const(ctg.clone());
            assert_eq!(dummy.is_record_ok(&mut read), result)
        }
    }
}
