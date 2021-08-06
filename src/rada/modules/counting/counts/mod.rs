use rust_htslib::bam::record::Record;

pub use locus::LocusCounts;
pub use strandedbuf::StrandedCountsBuffer;
pub use unstrandedbuf::UnstrandedCountsBuffer;

mod locus;
mod strandedbuf;
mod unstrandedbuf;

pub struct CountsBufferContent<'a> {
    pub forward: Option<&'a [LocusCounts]>,
    pub reverse: Option<&'a [LocusCounts]>,
    pub unstranded: Option<&'a [LocusCounts]>,
}

impl<'a> CountsBufferContent<'a> {
    pub fn capacity(&self) -> usize {
        self.forward.map_or(0, |x| x.len())
            + self.reverse.map_or(0, |x| x.len())
            + self.unstranded.map_or(0, |x| x.len())
    }
}

pub trait CountsBuffer {
    fn resize(&mut self, len: u32) -> &mut Self;
    // record must NOT be mutable here, yet some const(!) methods in rust_htslib require mutable(!) instance
    fn buffer_for(&mut self, record: &mut Record) -> &mut [LocusCounts];
    fn counts(&self) -> CountsBufferContent;
    fn len(&self) -> u32;
    fn reset(&mut self) -> &mut Self;
}
