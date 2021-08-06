use rust_htslib::bam::record::Record;

use super::{CountsBuffer, LocusCounts};
use crate::rada::modules::counting::counts::CountsBufferContent;

pub struct UnstrandedCountsBuffer {
    buffer: Vec<LocusCounts>,
}

impl UnstrandedCountsBuffer {
    pub(crate) fn new(maxsize: u32) -> Self {
        let mut buffer = Vec::new();
        buffer.reserve(maxsize as usize);
        UnstrandedCountsBuffer { buffer }
    }
}

impl CountsBuffer for UnstrandedCountsBuffer {
    fn resize(&mut self, len: u32) -> &mut Self {
        let len = len as usize;
        if self.buffer.len() != len {
            self.buffer.resize(len, LocusCounts::zeros());
        }
        self
    }

    #[inline]
    fn buffer_for(&mut self, _: &mut Record) -> &mut [LocusCounts] {
        &mut self.buffer
    }

    fn counts(&self) -> CountsBufferContent {
        CountsBufferContent {
            forward: None,
            reverse: None,
            unstranded: Some(&self.buffer),
        }
    }

    fn len(&self) -> u32 {
        self.buffer.len() as u32
    }

    fn reset(&mut self) -> &mut Self {
        self.buffer.fill(LocusCounts::zeros());
        self
    }
}
