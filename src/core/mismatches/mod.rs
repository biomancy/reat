use std::ops::Range;

use bio_types::genome::AbstractInterval;
use bio_types::strand::Strand;

use crate::core::rpileup::ncounters::AggregatedNucCounts;

pub mod prefilters;
pub mod roi;
pub mod site;

pub struct Restranded<T> {
    pub unknown: Option<T>,
    pub forward: Option<T>,
    pub reverse: Option<T>,
}

pub trait BatchedMismatches: AbstractInterval + Sized {
    type Flattened;

    fn strand(&self) -> Strand;

    // Mutability
    fn filter(self, mask: Vec<bool>) -> Self;
    fn restrand(self, strands: Vec<Strand>) -> Restranded<Self>;
    fn restrand_all(&mut self, strand: Strand);

    // Once batched representation is not needed, results might be flattened
    // (object of arrays -> array of objects)
    fn flatten(self) -> Vec<Self::Flattened>;
}

#[derive(Clone)]
pub struct MismatchesIntermediate<T: BatchedMismatches> {
    pub retained: Vec<T>,
    pub other: Vec<T>,
}

pub trait BatchedMismatchesBuilder<'a, NucCounts: AggregatedNucCounts<'a>> {
    type Mismatches: BatchedMismatches;

    fn build(&mut self, nc: NucCounts) -> MismatchesIntermediate<Self::Mismatches>;
}
