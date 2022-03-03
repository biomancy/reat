use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;

pub use batched::REATBatchedROIMismatches;
pub use builder::REATROIMismatchesBuilder;
pub use flat::REATROIMismatches;
pub use msummary::MismatchesSummary;

use crate::core::dna::NucCounts;
use crate::core::mismatches::BatchedMismatches;
use crate::core::workload::ROI;

mod batched;
mod builder;
mod flat;
mod msummary;

pub type ROIMismatchesPreview = MismatchesSummary;

pub trait ROIMismatches {
    // Original ROI record
    fn roi(&self) -> &ROI;
    // Transcription strand
    fn trstrand(&self) -> Strand;
    // Number of unique fragments covering the ROI
    fn coverage(&self) -> u32;
    // Total nucleotides in the given ROI (after masking)
    fn prednuc(&self) -> NucCounts;
    // Observed mismatches relative to the predicted reference
    fn mismatches(&self) -> &MismatchesSummary;
}

pub trait BatchedROIMismatches: BatchedMismatches
where
    Self::Flattened: ROIMismatches,
{
    // Original bed intervals
    fn premasked(&self) -> &[Range<Position>];
    // Minimal enclosing interval for each ROI after all masking applied
    fn postmasked(&self) -> &[Range<Position>];
    // Sub-intervals for each ROI after all masking applied
    fn subintervals(&self) -> &[Vec<Range<Position>>];
    // Number of unique reads covering the given region
    fn coverage(&self) -> &[u32];
    // Number of predicted nucleotides in each ROI;
    fn prednuc(&self) -> &[NucCounts];
    // Number of mismatches between predicted and sequenced nucleotides
    fn mismatches(&self) -> &[MismatchesSummary];
}
