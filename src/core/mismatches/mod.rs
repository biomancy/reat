use bio_types::genome::AbstractInterval;
use bio_types::strand::Strand;

use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::strandutil::StrandedData;

pub mod prefilters;
pub mod roi;
pub mod site;
pub mod utils;

pub type StrandingCounts = StrandedData<usize>;

pub trait BatchedMismatches: AbstractInterval + Sized {
    type Flattened;

    fn trstrand(&self) -> Strand;

    // Mutability
    fn filter(self, mask: Vec<bool>) -> Self;
    fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>>;
    fn restrand_all(&mut self, strand: Strand);

    // Once batched representation is not needed, results might be flattened
    // (object of arrays -> array of objects)
    fn flatten(self) -> Vec<Self::Flattened>;
}

#[derive(Clone)]
pub struct FilteredBatchedMismatches<T: BatchedMismatches> {
    pub retained: Vec<T>,
    pub other: Vec<T>,
}

pub trait BatchedMismatchesBuilder<'a, NucCounts: AggregatedNucCounts<'a>> {
    type Mismatches: BatchedMismatches;

    fn build(&mut self, nc: NucCounts) -> FilteredBatchedMismatches<Self::Mismatches>;
}
