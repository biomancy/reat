use bio_types::genome::Interval;
#[cfg(test)]
use mockall::{automock, mock};

pub use basecnt::BaseNucCounter;
pub use strandcnt::StrandedNucCounter;

use crate::core::counting::buffers::{CountsBuffer, NucCounts};
use crate::core::counting::buffers::{IntervalCounts, RawCounts};
use crate::core::read::AlignedRead;
use crate::core::workload::ROIWorkload;
use bio_types::strand::Strand;

mod basecnt;
mod strandcnt;

#[derive(Debug, PartialEq, Clone)]
pub struct CountingResult<'a> {
    pub name: &'a str,
    pub strand: Strand,
    pub roi: &'a Interval,
    pub cnts: RawCounts<'a>,
}

#[cfg_attr(test, automock)]
pub trait NucCounter<R: AlignedRead> {
    // Buffer details
    fn interval(&self) -> &Interval;
    fn rois(&self) -> &[Interval];
    // Number of counted nucleotides for the forward strand
    fn counted<'a>(&'a mut self) -> &[NucCounts];
    // The total number of mapped reads that passed all filters and were counted in at least 1 position
    fn mapped(&self) -> u32;
    fn empty(&self) -> bool {
        self.mapped() == 0
    }
    // Count read R
    fn count(&mut self, read: &mut R);
    // Finalize counting and return results. May be called more than once(i.e. iteratively)
    fn results<'a>(&'a self) -> Vec<CountingResult<'a>>;
    // Reset buffers
    fn reset(&mut self, workload: ROIWorkload);
}
