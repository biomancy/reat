use bio_types::strand::ReqStrand;
use rust_htslib::bam::record::Record;

use crate::rada::modules::stranding::deduct::StrandDeductor;

use super::{CountsBuffer, CountsBufferContent, LocusCounts};

pub struct StrandedCountsBuffer<Deductor: StrandDeductor> {
    forward: Vec<LocusCounts>,
    reverse: Vec<LocusCounts>,
    strand_deductor: Deductor,
}

impl<Deductor: StrandDeductor> StrandedCountsBuffer<Deductor> {
    fn new(maxsize: u32, strand_deductor: Deductor) -> Self {
        let (mut plstrand, mut mnstrand) = (Vec::new(), Vec::new());
        plstrand.reserve(maxsize as usize);
        mnstrand.reserve(maxsize as usize);
        StrandedCountsBuffer {
            forward: plstrand,
            reverse: mnstrand,
            strand_deductor,
        }
    }
}

impl<Deductor: StrandDeductor> CountsBuffer for StrandedCountsBuffer<Deductor> {
    fn resize(&mut self, len: u32) -> &mut Self {
        let maxsize = len as usize;
        if self.forward.len() != maxsize {
            self.forward.resize(maxsize, LocusCounts::zeros());
            self.reverse.resize(maxsize, LocusCounts::zeros());
        }
        self
    }

    #[inline]
    fn buffer_for(&mut self, record: &mut Record) -> &mut [LocusCounts] {
        match self.strand_deductor.deduce(record) {
            ReqStrand::Forward => &mut self.forward,
            ReqStrand::Reverse => &mut self.reverse,
        }
    }

    fn counts(&self) -> CountsBufferContent {
        CountsBufferContent {
            forward: Some(&self.forward),
            reverse: Some(&self.reverse),
            unstranded: None,
        }
    }

    fn len(&self) -> u32 {
        self.forward.len() as u32
    }

    fn reset(&mut self) -> &mut Self {
        self.forward.fill(LocusCounts::zeros());
        self.reverse.fill(LocusCounts::zeros());
        self
    }
}
