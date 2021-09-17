use std::path::PathBuf;

use bio_types::genome::Interval;
use bio_types::strand::Strand;
#[cfg(test)]
use mockall::automock;

pub use base::BaseRunner;

use crate::core::counting::buffers::{CountsBuffer, RawCounts};
use crate::core::dna::Nucleotide;
use crate::core::filtering::summary::{LocusSummaryFilter, ROISummaryFilter};
use crate::core::stranding::predict::{LocusStrandPredictor, ROIStrandPredictor};
use crate::core::summary::{LocusSummary, ROISummary};
use crate::core::workload::ROIWorkload;

mod base;

pub struct RunResult<'a> {
    pub interval: &'a Interval,
    pub name: &'a str,
    pub strand: Strand,
    pub cnts: RawCounts<'a>,
    pub reference: &'a [Nucleotide],
}

#[cfg_attr(test, automock)]
pub trait Runner {
    fn run<'a>(&'a mut self, workload: ROIWorkload) -> Vec<RunResult<'a>>;
    fn mapped(&self) -> u32;
}
