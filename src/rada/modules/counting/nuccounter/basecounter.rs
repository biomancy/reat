use std::cmp::min;

use bio_types::genome::{AbstractInterval, Interval};
use bio_types::sequence::SequenceRead;
use derive_more::Constructor;
use rust_htslib::bam::record::{Cigar, Record};

use crate::rada::modules::filtering::reads::ReadsFilter;

use super::super::counts::{CountsBuffer, CountsBufferContent};
use super::{NucCounter, NucCounterContent};
use crate::rada::modules::read::AlignedRead;
use std::marker::PhantomData;

pub struct BaseNucCounter<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>> {
    filter: Filter,
    buffer: Buffer,
    roi: Interval,
    phantom: PhantomData<R>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>>
    BaseNucCounter<R, Filter, Buffer>
{
    fn is_record_ok(&self, record: &mut R) -> bool {
        !self.filter.is_read_ok(record) || record.contig() != self.roi.contig()
    }

    pub fn new(filter: Filter, mut buffer: Buffer, roi: Interval) -> Self {
        buffer.reset((roi.range().end - roi.range().start) as u32);
        BaseNucCounter {
            filter,
            buffer,
            roi,
            phantom: Default::default(),
        }
    }
}

impl<R: AlignedRead, Filter: ReadsFilter<R>, Buffer: CountsBuffer<R>> NucCounter<R>
    for BaseNucCounter<R, Filter, Buffer>
{
    fn roi(&self) -> &Interval {
        &self.roi
    }

    // TODO: Refactor and unit-test this.
    fn process(&mut self, record: &mut R) {
        if !self.is_record_ok(record) {
            return;
        }
        let counts = self.buffer.buffer_for(record);

        let sequence = record.seq();
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
                        } else if self.filter.is_base_ok(record, seqpos) {
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

    fn reset(&mut self, roi: Interval) {
        let size = (roi.range().end - roi.range().start) as u32;
        self.roi = roi;
        self.buffer.reset(size);
    }

    fn content<'a>(&'a self) -> NucCounterContent<'a> {
        NucCounterContent {
            interval: self.roi.clone(),
            counts: self.buffer.content(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use crate::rada::modules::filtering::reads::MockReadsFilter;

    use super::super::super::counts::MockCountsBuffer;
    use super::*;

    //
    // #[cfg(test)]
    // fn reset() {
    //     let mut counter = BaseNucCounter::<>::new(
    //         MockReadsFilter::new(),
    //         MockCountsBuffer::new(),
    //         Interval::new("".into(), 1..2),
    //     );
    //     assert_eq!(counter.buffer.len(), 1);
    //
    //     let newroi = Interval::new("".to_string(), 100..200);
    //     counter.reset(newroi);
    //     assert_eq!(counter.buffer.len(), 100);
    // }
    //
    // #[cfg(test)]
    // fn content() {
    //     let mut buffer = MockCountsBuffer::new();
    //     buffer
    //         .expect_content()
    //         .once()
    //         .return_const(CountsBufferContent {
    //             forward: None,
    //             reverse: None,
    //             unstranded: None,
    //         });
    //     let counter = BaseNucCounter::new(
    //         MockReadsFilter::new(),
    //         buffer,
    //         Interval::new("".into(), 1..2),
    //     );
    //
    //     let content = counter.content();
    //     assert_eq!(content.interval, counter.roi);
    //     assert_eq!(content.counts.forward, None);
    //     assert_eq!(content.counts.reverse, None);
    //     assert_eq!(content.counts.unstranded, None);
    // }
}
