use std::cmp::min;
use std::marker::PhantomData;
use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{ReqStrand, Strand};
use rust_htslib::bam::record::Cigar;

use crate::core::dna::NucCounts;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;

#[derive(Clone)]
pub struct BaseNucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    // Filtering parameters
    trim5: usize,
    trim3: usize,
    rfilter: Filter,
    // Caches
    buffer: Vec<NucCounts>,
    matched: Vec<Range<u32>>,
    mapped: u32,
    // Current interval
    interval: Interval,
    phantom: PhantomData<fn() -> R>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>> BaseNucCounter<R, Filter> {
    pub fn new(maxbuf: usize, filter: Filter, trim5: u16, trim3: u16) -> Self {
        BaseNucCounter {
            rfilter: filter,
            interval: Interval::new("".to_string(), 0..0),
            buffer: Vec::with_capacity(maxbuf),
            matched: Vec::with_capacity(20),
            mapped: 0,
            trim5: trim5 as usize,
            trim3: trim3 as usize,
            phantom: Default::default(),
        }
    }

    #[inline]
    pub fn interval(&self) -> &Interval {
        &self.interval
    }

    #[inline]
    pub fn counted(&self) -> &[NucCounts] {
        &self.buffer
    }

    #[inline]
    pub fn mapped(&self) -> u32 {
        self.mapped
    }

    #[inline]
    pub fn reset(&mut self, interval: Interval) {
        let newlen = interval.range().end - interval.range().start;
        debug_assert!(newlen > 0);
        self.buffer.clear();
        self.buffer.resize(newlen as usize, NucCounts::zeros());

        self.mapped = 0;
        self.interval = interval;
    }

    pub fn count(&mut self, read: &R) -> &[Range<u32>] {
        self.matched.clear();

        if self.is_record_ok(read) {
            self.implprocess(read);

            if !self.matched.is_empty() {
                self.mapped += 1;
            }
        }
        return &self.matched;
    }

    #[inline]
    fn is_record_ok(&self, record: &R) -> bool {
        self.rfilter.is_read_ok(record) && record.contig() == self.interval.contig()
    }

    fn implprocess(&mut self, read: &R) {
        let sequence = read.seq();

        let (mut roipos, mut seqpos) = (read.pos() - self.interval.range().start as i64, 0usize);
        let roisize = (self.interval.range().end - self.interval.range().start) as i64;

        // Read is too short
        if read.len() <= (self.trim5 + self.trim3) {
            return;
        }

        let (minseqpos, maxseqpos) = match read.strand() {
            ReqStrand::Forward => (self.trim5, read.len() - self.trim3),
            ReqStrand::Reverse => (self.trim3, read.len() - self.trim5),
        };

        for block in read.cigar().iter() {
            if roipos >= roisize || seqpos >= maxseqpos {
                break;
            }
            match block {
                Cigar::Match(ops) | Cigar::Equal(ops) | Cigar::Diff(ops) => {
                    // fast-end when possible
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
                        if seqpos >= minseqpos && seqpos < maxseqpos && self.rfilter.is_base_ok(read, seqpos) {
                            debug_assert!(roipos >= 0);
                            // From the SAM specification: No assumptions can be made on the letter cases
                            match sequence[seqpos as usize] {
                                b'A' | b'a' => self.buffer[roipos as usize].A += 1,
                                b'T' | b't' => self.buffer[roipos as usize].T += 1,
                                b'G' | b'g' => self.buffer[roipos as usize].G += 1,
                                b'C' | b'c' => self.buffer[roipos as usize].C += 1,
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
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use rust_htslib::bam::record::CigarString;

    use shortcats::*;

    use crate::core::dna::NucCounts;
    use crate::core::read::MockRead;
    use crate::core::rpileup::ncounters::filters::MockReadsFilter;

    use super::*;

    mod shortcats {
        #![allow(non_snake_case)]

        use rust_htslib::bam::record::Cigar;
        use rust_htslib::bam::record::Cigar::*;

        use crate::core::dna::NucCounts;

        pub fn Z() -> NucCounts {
            NucCounts::zeros()
        }
        pub fn A() -> NucCounts {
            NucCounts::A(1)
        }
        pub fn C() -> NucCounts {
            NucCounts::C(1)
        }
        pub fn G() -> NucCounts {
            NucCounts::G(1)
        }
        pub fn T() -> NucCounts {
            NucCounts::T(1)
        }

        pub fn M(x: u32) -> Cigar {
            Match(x)
        }
        pub fn X(x: u32) -> Cigar {
            Diff(x)
        }
        pub fn E(x: u32) -> Cigar {
            Equal(x)
        }
        pub fn D(x: u32) -> Cigar {
            Del(x)
        }
        pub fn N(x: u32) -> Cigar {
            RefSkip(x)
        }
        pub fn H(x: u32) -> Cigar {
            HardClip(x)
        }
        pub fn P(x: u32) -> Cigar {
            Pad(x)
        }
        pub fn S(x: u32) -> Cigar {
            SoftClip(x)
        }
        pub fn I(x: u32) -> Cigar {
            Ins(x)
        }
    }

    fn run(
        trim: (u16, u16),
        roi: Range<u64>,
        pos: i64,
        seq: &str,
        strand: ReqStrand,
        okbase: Vec<bool>,
        cigar: Vec<Cigar>,
        excounts: &[NucCounts],
        exmatch: &[Range<u32>],
    ) {
        let roi = Interval::new("".into(), roi);
        let roisize = roi.range().end - roi.range().start;

        let mut filter = MockReadsFilter::new();
        for isok in okbase {
            filter.expect_is_base_ok().once().return_const(isok);
        }

        let mut counter = BaseNucCounter::new((roisize + 1) as usize, filter, trim.0, trim.1);
        counter.reset(roi);

        let mut read = MockRead::new();
        read.expect_pos().return_const(pos);
        read.expect_len().return_const(seq.len());
        read.expect_cigar().return_once(move || CigarString(cigar).into_view(pos));
        read.expect_strand().return_const(strand);
        let seq = String::from(seq);
        read.expect_seq().returning(move || seq.as_bytes().to_vec());

        counter.implprocess(&mut read);

        assert_eq!(counter.buffer, excounts);
        assert_eq!(counter.matched, exmatch);
    }

    #[test]
    #[allow(non_snake_case)]
    fn implprocess() {
        // Query consuming operations only
        for op in [M, X, E] {
            // complete overlap with the region
            run(
                (0, 0),
                0..4,
                0,
                "ACGT",
                ReqStrand::Forward,
                vec![true, true, true, true],
                vec![op(4)],
                &vec![A(), C(), G(), T()],
                &vec![0..4],
            );
            // inside region
            run(
                (0, 0),
                0..4,
                1,
                "AC",
                ReqStrand::Forward,
                vec![true, true],
                vec![op(2)],
                &vec![Z(), A(), C(), Z()],
                &vec![1..3],
            );
            // completely out of the region
            run((0, 0), 0..4, 4, "ACGT", ReqStrand::Reverse, vec![], vec![op(4)], &vec![Z(), Z(), Z(), Z()], &vec![]);
            // end out of the region
            run(
                (0, 0),
                0..4,
                2,
                "ACGT",
                ReqStrand::Forward,
                vec![true, true],
                vec![op(4)],
                &vec![Z(), Z(), A(), C()],
                &vec![2..4],
            );
            // start out of the region
            run(
                (0, 0),
                0..4,
                -1,
                "ACGT",
                ReqStrand::Reverse,
                vec![true, true, true],
                vec![op(4)],
                &vec![C(), G(), T(), Z()],
                &vec![0..3],
            );
        }

        // No-ops + reference/query consuming operations only
        let empty = vec![Z()].repeat(4);
        for op in [D, N, H, P, S, I] {
            // complete overlap with the region
            run((0, 0), 0..4, 0, "ACGT", ReqStrand::Forward, vec![], vec![op(4)], &empty, &vec![]);
            // inside region
            run((0, 0), 0..4, 1, "AC", ReqStrand::Reverse, vec![], vec![op(2)], &empty, &vec![]);
            // completely out of the region
            run((0, 0), 0..4, 4, "ACGT", ReqStrand::Forward, vec![], vec![op(4)], &empty, &vec![]);
            // end out of the region
            run((0, 0), 0..4, 2, "ACGT", ReqStrand::Reverse, vec![], vec![op(4)], &empty, &vec![]);
            // start out of the region
            run((0, 0), 0..4, -1, "ACGT", ReqStrand::Forward, vec![], vec![op(4)], &empty, &vec![]);
        }

        // Complex queries
        run(
            (0, 0),
            2..5,
            0,
            "AGC",
            ReqStrand::Forward,
            vec![true],
            vec![D(4), M(1), N(3), M(2)],
            &vec![Z(), Z(), A()],
            &vec![2..3],
        );
        run((0, 0), 2..5, 0, "AGC", ReqStrand::Reverse, vec![false], vec![M(3)], &vec![Z(), Z(), Z()], &vec![]);
        run(
            (0, 0),
            2..5,
            0,
            "AGC",
            ReqStrand::Reverse,
            vec![true, true],
            vec![M(1), N(2), M(2)],
            &vec![Z(), G(), C()],
            &vec![1..3],
        );
        run(
            (0, 0),
            1..5,
            0,
            "NNN",
            ReqStrand::Forward,
            vec![true, true],
            vec![M(1), N(2), M(2)],
            &vec![Z(), Z(), Z(), Z()],
            &vec![2..4],
        );

        run(
            (0, 0),
            2..5,
            0,
            "AGC",
            ReqStrand::Forward,
            vec![true, true],
            vec![D(2), M(1), N(1), M(1), N(1), M(1)],
            &vec![A(), Z(), G()],
            &vec![0..1, 2..3],
        );
        run(
            (0, 0),
            0..10,
            0,
            "ACGTNNGTAG",
            ReqStrand::Reverse,
            vec![true, true, false, true, true, true, true, false, true, true],
            vec![M(10)],
            &vec![A(), C(), Z(), T(), Z(), Z(), G(), Z(), A(), G()],
            &vec![0..2, 3..7, 8..10],
        );
    }

    #[test]
    fn trim() {
        // normal trimming
        run(
            (1, 2),
            1..6,
            1,
            "ACGT",
            ReqStrand::Forward,
            vec![true],
            vec![M(4)],
            &vec![Z(), C(), Z(), Z(), Z()],
            &vec![1..2],
        );
        run(
            (1, 2),
            1..6,
            1,
            "ACGT",
            ReqStrand::Reverse,
            vec![true],
            vec![M(4)],
            &vec![Z(), Z(), G(), Z(), Z()],
            &vec![2..3],
        );

        // trim all bases from 5` or 3`
        let ex = (vec![Z()].repeat(10), vec![]);
        for trim in [3, 6] {
            run((trim, 0), 10..20, 10, "ACG", ReqStrand::Forward, vec![], vec![M(3)], &ex.0, &ex.1);
            run((0, trim), 10..20, 17, "ACG", ReqStrand::Forward, vec![], vec![M(3)], &ex.0, &ex.1);
            run((trim, trim), 10..20, 17, "ACG", ReqStrand::Reverse, vec![], vec![M(3)], &ex.0, &ex.1);
        }

        // Read out of the roi + trim
        run((2, 1), 3..4, 0, "ACGTA", ReqStrand::Forward, vec![true], vec![M(5)], &vec![T()], &vec![0..1]);
        run((2, 1), 3..4, 1, "ACGTA", ReqStrand::Reverse, vec![true], vec![M(5)], &vec![G()], &vec![0..1]);

        run((0, 3), 0..2, 0, "AGTA", ReqStrand::Forward, vec![true], vec![M(4)], &vec![A(), Z()], &vec![0..1]);
        run((0, 3), 0..2, 0, "AGTA", ReqStrand::Reverse, vec![], vec![M(4)], &vec![Z(), Z()], &vec![]);

        run((2, 0), 2..4, 1, "CGTA", ReqStrand::Forward, vec![true], vec![M(4)], &vec![Z(), T()], &vec![1..2]);
        run((2, 0), 2..4, 1, "CGTA", ReqStrand::Reverse, vec![true], vec![M(4)], &vec![G(), Z()], &vec![0..1]);
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
            let mut dummy = BaseNucCounter::new(1, filter, 4, 0);
            dummy.reset(Interval::new(contig.clone(), 0..1));

            let mut read = MockRead::new();
            read.expect_contig().return_const(ctg.clone());
            assert_eq!(dummy.is_record_ok(&mut read), result)
        }
    }
}
