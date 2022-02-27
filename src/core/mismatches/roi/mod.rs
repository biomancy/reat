use std::borrow::Borrow;
use std::ops::Range;
use std::path::Path;

use bio_types::genome::{Interval, Position};
use bio_types::strand::Strand;

pub use batched::REATBinnedROIMismatches;
pub use flat::REATROIMismatches;
pub use msummary::MismatchesSummary;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::mismatches::BatchedMismatches;
use crate::core::workload::ROI;

mod batched;
mod flat;
mod msummary;

pub type ROIMismatchesPreview = MismatchesSummary;

pub trait ROIMismatches {
    fn roi(&self) -> &ROI;
    fn strand(&self) -> Strand;
    fn coverage(&self) -> u32;
    fn prednuc(&self) -> NucCounts;
    fn mismatches(&self) -> &MismatchesSummary;
}

pub trait BinnedROIMismatches: BatchedMismatches
where
    Self::Flattened: ROIMismatches,
{
    // Minimal enclosing interval for each roi
    fn rois(&self) -> &[Range<Position>];
    // Rois' pieces after filtering excluded regions
    fn pieces(&self) -> &[Vec<Range<Position>>];
    // Original bed intervals
    fn premasked(&self) -> &[Range<Position>];
    // Number of unique reads covering the given region
    fn coverage(&self) -> &[u32];
    // Number of predicted nucleotides in each ROI;
    fn prednuc(&self) -> &[NucCounts];
    // Number of mismatches between predicted and sequenced nucleotides
    fn mismatches(&self) -> &[MismatchesSummary];
}
