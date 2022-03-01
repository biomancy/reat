use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{ReqStrand, Strand};

pub use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;

pub mod cnt;
pub mod filters;

pub trait NucCounter<'a, R: AlignedRead>: ReadsCollider<'a, R>
where
    Self::ColliderResult: AggregatedNucCounts<'a>,
{
}

pub trait AggregatedNucCounts<'a>: AbstractInterval {
    type ItemInfo;

    fn items(&'a self) -> &'a [AggregatedNucCountsItem<'a, Self::ItemInfo>];
    fn items_mut<'b>(&'b mut self) -> &'b mut [AggregatedNucCountsItem<'a, Self::ItemInfo>];
    fn consume(self) -> (&'a str, Vec<AggregatedNucCountsItem<'a, Self::ItemInfo>>);
}

pub struct AggregatedNucCountsItem<'a, ItemInfo> {
    pub info: ItemInfo,
    pub range: Range<Position>,
    pub forward: Option<&'a [NucCounts]>,
    pub reverse: Option<&'a [NucCounts]>,
    pub unstranded: Option<&'a [NucCounts]>,
}

impl<'a, ItemInfo> AggregatedNucCountsItem<'a, ItemInfo> {
    pub fn seqnuc<'b>(&'a self, buffer: &'b mut Vec<NucCounts>) -> Option<&'a [NucCounts]> {
        // Gather sequenced nucleotides in each position
        match (self.forward, self.reverse, self.unstranded) {
            (Some(c), None, None) => Some(c),
            (None, Some(c), None) => Some(c),
            (None, None, Some(c)) => Some(c),
            (_, _, _) => {
                buffer.clear();
                buffer.resize(self.range.end as usize - self.range.start as usize, Default::default());
                for x in [self.forward, self.reverse, self.unstranded] {
                    if let Some(x) = x {
                        debug_assert!(x.len() == buffer.len());
                        for (bc, xc) in buffer.iter_mut().zip(x.iter()) {
                            *bc += *xc;
                        }
                    }
                }
                None
            }
        }
    }
}
