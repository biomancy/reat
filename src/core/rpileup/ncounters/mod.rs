use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{ReqStrand, Strand};

pub use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::mismatches::BatchedMismatches;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;

pub mod cnt;
pub mod filters;

pub trait NucCounter<'a, R: AlignedRead>: ReadsCollider<'a, R>
where
    Self::ColliderResult: AggregatedNucCounts<'a>,
{
}

pub struct AggregatedNucCountsItem<'a, Meta> {
    pub meta: Meta,
    pub range: Range<Position>,
    pub forward: Option<&'a [NucCounts]>,
    pub reverse: Option<&'a [NucCounts]>,
    pub unstranded: Option<&'a [NucCounts]>,
}

pub trait AggregatedNucCounts<'a>: AbstractInterval {
    type Meta;

    fn items(&'a self) -> &'a [AggregatedNucCountsItem<'a, Self::Meta>];
    fn consume(self) -> Vec<AggregatedNucCountsItem<'a, Self::Meta>>;
}
