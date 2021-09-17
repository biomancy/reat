use std::cmp::{max, min};
use std::marker::PhantomData;
use std::ops::{Index, Range};

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::Strand;
use itertools::Itertools;
use rust_htslib::bam::record::Cigar;

use crate::core::counting::buffers::IntervalCounts;
use crate::core::counting::buffers::NucCounts;
use crate::core::counting::nuccounter::{CountingResult, NucCounter};
use crate::core::filtering::reads::ReadsFilter;
use crate::core::read::AlignedRead;

use super::super::buffers::CountsBuffer;
use crate::core::workload::ROIWorkload;

#[derive(Clone)]
pub struct BaseNucCounter<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer> {
    filter: Filter,
    buffer: Buffer,
    matched: Vec<Range<u32>>,
    mapped: u32,
    phantom: PhantomData<fn() -> R>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer> BaseNucCounter<R, Filter, Buffer> {
    #[inline]
    fn is_record_ok(&self, record: &mut R) -> bool {
        self.filter.is_read_ok(record) && record.contig() == self.buffer.interval().contig()
    }

    pub fn new(filter: Filter, buffer: Buffer) -> Self {
        let matched = Vec::with_capacity(5);
        BaseNucCounter { filter, buffer, matched, mapped: 0, phantom: Default::default() }
    }

    fn implprocess(&mut self, read: &mut R) -> &[NucCounts] {
        let sequence = read.seq();

        let (mut roipos, mut seqpos) = (read.pos() - self.buffer.interval().range().start as i64, 0usize);
        let roisize = (self.buffer.interval().range().end - self.buffer.interval().range().start) as i64;

        let counts = self.buffer.buffer_mut();
        self.matched.clear();

        for cigar in read.cigar().iter() {
            if roipos >= roisize {
                break;
            }
            match cigar {
                Cigar::Match(ops) | Cigar::Equal(ops) | Cigar::Diff(ops) => {
                    let end = min(roisize as i64, roipos + *ops as i64);
                    // fast-forward when possible
                    if roipos < 0 {
                        let skip = min(roipos.abs() as u32, *ops);
                        roipos += skip as i64;
                        seqpos += skip as usize;
                    }

                    let mut prevmatched: Option<u32> = None;
                    let start = roipos;
                    for _ in start..end {
                        debug_assert!(roipos < roisize);
                        if self.filter.is_base_ok(read, seqpos) {
                            debug_assert!(roipos >= 0);
                            // From the SAM specification: No assumptions can be made on the letter cases
                            match sequence[seqpos as usize] {
                                b'A' | b'a' => counts[roipos as usize].A += 1,
                                b'T' | b't' => counts[roipos as usize].T += 1,
                                b'G' | b'g' => counts[roipos as usize].G += 1,
                                b'C' | b'c' => counts[roipos as usize].C += 1,
                                _ => {}
                            }
                            if prevmatched.is_none() {
                                prevmatched = Some(roipos as u32);
                            }
                        } else if let Some(m) = prevmatched {
                            self.matched.push(m..roipos as u32);
                            prevmatched = None;
                        }
                        roipos += 1;
                        seqpos += 1;
                    }
                    if let Some(m) = prevmatched {
                        self.matched.push(m..roipos as u32);
                        prevmatched = None;
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

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer> NucCounter<R> for BaseNucCounter<R, Filter, Buffer> {
    #[inline]
    fn interval(&self) -> &Interval {
        self.buffer.interval()
    }

    #[inline]
    fn rois(&self) -> &[Interval] {
        self.buffer.rois()
    }

    #[inline]
    fn counted(&mut self) -> &[NucCounts] {
        self.buffer.buffer()
    }

    #[inline]
    fn mapped(&self) -> u32 {
        self.mapped
    }

    fn count(&mut self, read: &mut R) {
        if !self.is_record_ok(read) {
            return;
        }
        self.implprocess(read);
        if !self.matched.is_empty() {
            self.buffer.add_matched(&self.matched);
            self.mapped += 1;
        }
    }

    fn results(&self) -> Vec<CountingResult> {
        self.buffer
            .results()
            .into_iter()
            .map(|x| CountingResult { name: x.name, strand: Strand::Unknown, roi: x.roi, cnts: x.cnts })
            .collect()
    }

    #[inline]
    fn reset(&mut self, info: ROIWorkload) {
        self.buffer.reset(info);
        self.mapped = 0;
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use mockall::predicate;
    use rust_htslib::bam::record::Cigar::*;
    use rust_htslib::bam::record::CigarString;

    use crate::core::counting::NucCounts;
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
        excounts: &[NucCounts],
        exmatch: &[Range<u32>],
    ) {
        let roi = Interval::new("".into(), roi);
        let roisize = roi.range().end - roi.range().start;

        let mut buffer = MockCountsBuffer::default();
        buffer.expect_reset().once().return_const(());
        buffer.expect_interval().return_const(roi.clone());
        buffer.expect_buffer_mut().once().returning(move || vec![NucCounts::zeros()].repeat(roisize as usize));

        let mut filter = MockReadsFilter::new();
        for isok in okbase {
            filter.expect_is_base_ok().once().return_const(isok);
        }

        let mut counter = BaseNucCounter::new(filter, buffer);
        counter.reset(ROIWorkload::default());

        let mut read = MockRead::new();
        read.expect_pos().return_const(pos);
        read.expect_cigar().return_once(move || CigarString(cigar).into_view(pos));
        let seq = String::from(seq);
        read.expect_seq().returning(move || seq.as_bytes().to_vec());

        assert_eq!(counter.implprocess(&mut read), excounts);
        assert_eq!(counter.matched, exmatch);
    }

    #[test]
    #[allow(non_snake_case)]
    fn implprocess() {
        let A = || NucCounts::A(1);
        let C = || NucCounts::C(1);
        let G = || NucCounts::G(1);
        let T = || NucCounts::T(1);
        let Z = || NucCounts::zeros();

        // Query consuming operations only
        let M = |x| Match(x);
        let X = |x| Diff(x);
        let E = |x| Equal(x);
        for op in [M, X, E] {
            // complete overlap with the region
            _implprocess(
                0..4,
                0,
                "ACGT",
                vec![true, true, true, true],
                vec![op(4)],
                &vec![A(), C(), G(), T()],
                &vec![0..4],
            );
            // inside region
            _implprocess(0..4, 1, "AC", vec![true, true], vec![op(2)], &vec![Z(), A(), C(), Z()], &vec![1..3]);
            // completely out of the region
            _implprocess(0..4, 4, "ACGT", vec![], vec![op(4)], &vec![Z(), Z(), Z(), Z()], &vec![]);
            // end out of the region
            _implprocess(0..4, 2, "ACGT", vec![true, true], vec![op(4)], &vec![Z(), Z(), A(), C()], &vec![2..4]);
            // start out of the region
            _implprocess(0..4, -1, "ACGT", vec![true, true, true], vec![op(4)], &vec![C(), G(), T(), Z()], &vec![0..3]);
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
            _implprocess(0..4, 0, "ACGT", vec![], vec![op(4)], &empty, &vec![]);
            // inside region
            _implprocess(0..4, 1, "AC", vec![], vec![op(2)], &empty, &vec![]);
            // completely out of the region
            _implprocess(0..4, 4, "ACGT", vec![], vec![op(4)], &empty, &vec![]);
            // end out of the region
            _implprocess(0..4, 2, "ACGT", vec![], vec![op(4)], &empty, &vec![]);
            // start out of the region
            _implprocess(0..4, -1, "ACGT", vec![], vec![op(4)], &empty, &vec![]);
        }

        // Complex queries
        _implprocess(
            2..5,
            0,
            "AGC",
            vec![true, true],
            vec![D(2), M(1), N(1), M(1), N(1), M(1)],
            &vec![A(), Z(), G()],
            &vec![0..1, 2..3],
        );
        _implprocess(2..5, 0, "AGC", vec![true], vec![D(4), M(1), N(3), M(2)], &vec![Z(), Z(), A()], &vec![2..3]);
        _implprocess(2..5, 0, "AGC", vec![false], vec![M(3)], &vec![Z(), Z(), Z()], &vec![]);
        _implprocess(2..5, 0, "AGC", vec![true, true], vec![M(1), N(2), M(2)], &vec![Z(), G(), C()], &vec![1..3]);
        _implprocess(1..5, 0, "NNN", vec![true, true], vec![M(1), N(2), M(2)], &vec![Z(), Z(), Z(), Z()], &vec![2..4]);

        _implprocess(
            0..10,
            0,
            "ACGTNNGTAG",
            vec![true, true, false, true, true, true, true, false, true, true],
            vec![M(10)],
            &vec![A(), C(), Z(), T(), Z(), Z(), G(), Z(), A(), G()],
            &vec![0..2, 3..7, 8..10],
        );
    }

    #[test]
    fn is_record_ok() {
        let contig = "".to_string();
        let interval = Interval::new(contig.clone(), 0..10);
        let wrong_contig = "!".to_string();
        for (isok, ctg, result) in [
            (true, &contig, true),
            (false, &contig, false),
            (true, &wrong_contig, false),
            (false, &wrong_contig, false),
        ] {
            let mut filter = MockReadsFilter::new();
            filter.expect_is_read_ok().once().return_const(isok);

            let mut buffer = MockCountsBuffer::default();
            buffer.expect_interval().return_const(interval.clone());
            let dummy = BaseNucCounter::new(filter, buffer);

            let mut read = MockRead::new();
            read.expect_contig().return_const(ctg.clone());
            assert_eq!(dummy.is_record_ok(&mut read), result)
        }
    }
}
