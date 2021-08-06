use std::cmp::min;
use std::fmt::Debug;

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::ReqStrand;
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::Read;

use crate::rada::modules::filtering::records::RecordsFilter;

use derive_more::Constructor;

use super::counts::{CountsBuffer, CountsBufferContent, LocusCounts};

pub struct NucCounterContent<'a> {
    pub interval: Interval,
    pub counts: CountsBufferContent<'a>,
}

pub trait NucCounter {
    fn process(&mut self, record: &mut Record);
    fn reset(&mut self, roi: &Interval);
    fn counts(&self) -> NucCounterContent;
}

#[derive(Constructor)]
pub struct BaseNucCounter<Filter: RecordsFilter, Buffer: CountsBuffer> {
    filter: Filter,
    buffer: Buffer,
    roi: Interval,
}

impl<Filter: RecordsFilter, Buffer: CountsBuffer> NucCounter for BaseNucCounter<Filter, Buffer> {
    fn process(&mut self, record: &mut Record) {
        if !self.filter.is_record_ok(record) || record.contig() != self.roi.contig() {
            return;
        }

        let counts = self.buffer.buffer_for(record);

        let sequence = record.seq().as_bytes();
        let qual = record.qual();

        let (mut offset, mut seqpos) = (record.pos() - self.roi.range().start as i64, 0usize);
        let size = self.roi.range().end - self.roi.range().start;
        for cigar in record.cigar().iter() {
            if offset >= size as i64 {
                break;
            }
            match cigar {
                Cigar::Match(ops) | Cigar::Equal(ops) | Cigar::Diff(ops) => {
                    // fast-forward when possible
                    let skip = if offset < 0 {
                        let skip = min(offset.abs() as u32, *ops);
                        offset += skip as i64;
                        seqpos += skip as usize;
                        skip
                    } else {
                        0
                    };
                    for _ in skip..*ops {
                        // out of bounds or position passes phread threshold
                        if offset >= size as i64 {
                            break;
                        } else if self.filter.is_base_ok(record, qual, seqpos) {
                            debug_assert!(offset >= 0);
                            match sequence[seqpos as usize] {
                                b'A' | b'a' => counts[offset as usize].A += 1,
                                b'T' | b't' => counts[offset as usize].T += 1,
                                b'G' | b'g' => counts[offset as usize].G += 1,
                                b'C' | b'c' => counts[offset as usize].C += 1,
                                base => {
                                    panic!("Unknown base {}", base)
                                }
                            }
                        }
                        offset += 1;
                        seqpos += 1;
                    }
                }
                Cigar::Del(ops) | Cigar::RefSkip(ops) => {
                    offset += *ops as i64;
                }
                Cigar::SoftClip(ops) | Cigar::Ins(ops) => {
                    seqpos += *ops as usize;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
    }

    fn reset(&mut self, roi: &Interval) {
        self.roi = roi.clone();
        self.buffer
            .resize((roi.range().end - roi.range().start) as u32)
            .reset();
    }

    fn counts(&self) -> NucCounterContent {
        NucCounterContent {
            interval: self.roi.clone(),
            counts: self.buffer.counts(),
        }
    }
}
