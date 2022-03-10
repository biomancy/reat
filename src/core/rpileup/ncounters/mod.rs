use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};

pub use crate::core::dna::NucCounts;
use crate::core::mismatches::StrandingCounts;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;
use crate::core::strandutil::StrandedData;

pub mod cnt;
pub mod cntitem;
pub mod filters;

pub trait NucCounter<'a, R: AlignedRead>: ReadsCollider<'a, R>
where
    Self::ColliderResult: AggregatedNucCounts<'a>,
{
}

pub trait AggregatedNucCounts<'a>: AbstractInterval {
    type ItemInfo;

    fn items(&self) -> &[AggregatedNucCountsItem<'a, Self::ItemInfo>];
    fn items_mut(&mut self) -> &mut [AggregatedNucCountsItem<'a, Self::ItemInfo>];
    fn consume(self) -> (&'a str, Vec<AggregatedNucCountsItem<'a, Self::ItemInfo>>);
}

pub struct AggregatedNucCountsItem<'a, ItemInfo> {
    pub info: ItemInfo,
    pub range: Range<Position>,
    pub mapped: StrandedData<u32>,
    pub counts: StrandedData<Option<&'a [NucCounts]>>,
}
