use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;

pub use builder::REATROIMismatchesBuilder;
pub use flat::REATROIMismatches;
pub use msummary::NucMismatches;
pub use vec::REATROIMismatchesVec;

use crate::core::dna::NucCounts;
use crate::core::mismatches::MismatchesVec;
use crate::core::workload::ROI;

mod builder;
mod flat;
mod msummary;
mod vec;

pub type ROIMismatchesPreview = NucMismatches;

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
    fn mismatches(&self) -> &NucMismatches;
}

pub trait ROIMismatchesVec: MismatchesVec
where
    Self::Flat: ROIMismatches,
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
    fn mismatches(&self) -> &[NucMismatches];
}

// use super::{AbstractInterval, StrandedData, StrandingCounts};
// #[cfg(test)]
// mock! {
//     pub BatchedROIMismatches {}
//     impl AbstractInterval for BatchedROIMismatches {
//         fn contig(&self) -> &str;
//         fn range(&self) -> Range<Position>;
//     }
//     impl BatchedMismatches for BatchedROIMismatches {
//         type Flattened = ();
//         fn trstrand(&self) -> Strand;
//         fn filter(self, mask: Vec<bool>) -> Self;
//         fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>>;
//         fn restrand_all(&mut self, strand: Strand);
//         fn flatten(self) -> Vec<<Self as BatchedMismatches>::Flattened>;
//     }
//     impl BatchedROIMismatches for BatchedROIMismatches {
//         fn premasked(&self) -> &[Range<Position>];
//         fn postmasked(&self) -> &[Range<Position>];
//         fn subintervals(&self) -> &[Vec<Range<Position>>];
//         fn coverage(&self) -> &[u32];
//         fn prednuc(&self) -> &[NucCounts];
//         fn mismatches(&self) -> &[MismatchesSummary];
//     }
// }
