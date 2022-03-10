use bio_types::genome::AbstractInterval;
use bio_types::strand::Strand;

use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::strandutil::StrandedData;

pub mod prefilters;
pub mod roi;
pub mod site;
pub mod utils;

pub type StrandingCounts = StrandedData<usize>;

pub trait MismatchesVec: AbstractInterval + Sized {
    type Flat;

    // Predicted transcription strand
    fn trstrand(&self) -> Strand;

    // Once grouped representation is not needed, results might be flattened
    // (object of arrays -> array of objects)
    fn flatten(self) -> Vec<Self::Flat>;
}

pub struct Context<T: MismatchesVec> {
    // Must be retained & printed no matter what
    pub retained: StrandedData<Option<T>>,
    // Other mismatches
    pub items: StrandedData<Option<T>>,
}

pub trait MismatchesBuilder<'a, NucCounts: AggregatedNucCounts<'a>> {
    type Mismatches: MismatchesVec;

    fn build(&'a mut self, nc: NucCounts) -> Context<Self::Mismatches>;
}
