use std::ops::Range;

use bio_types::genome::Interval;
#[cfg(test)]
use mockall::automock;

pub use crate::core::counting::nuccounts::NucCounts;
use crate::core::workload::ROIWorkload;

mod flatbuf;
mod roibuf;

pub use flatbuf::FlatBuffer;
pub use roibuf::ROIBuffer;

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct RawCounts<'a> {
    pub nuc: &'a [NucCounts],
    pub coverage: u32,
}

#[derive(Eq, PartialEq, Debug)]
pub struct IntervalCounts<'a> {
    pub roi: &'a Interval,
    pub name: &'a str,
    pub cnts: RawCounts<'a>,
}

#[cfg_attr(test, automock)]

pub trait CountsBuffer {
    // Content of the buffer
    fn interval(&self) -> &Interval;
    fn rois(&self) -> &[Interval];
    fn buffer(&self) -> &[NucCounts];
    // Counting methods
    fn buffer_mut(&mut self) -> &mut [NucCounts];
    fn add_matched(&mut self, blocks: &[Range<u32>]);
    // ROI results & reset
    #[allow(clippy::needless_lifetimes)]
    fn results<'a>(&'a self) -> Vec<IntervalCounts<'a>>;
    fn reset(&mut self, info: ROIWorkload);
}
