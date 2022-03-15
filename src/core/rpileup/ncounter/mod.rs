use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};

pub use crate::core::dna::NucCounts;
use crate::core::mismatches::StrandingCounts;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;
use crate::core::strandutil::Stranded;

pub mod cnt;
pub mod cntitem;
pub mod filters;

pub struct InnerNucCounts<'a, Data> {
    pub data: Data,
    pub range: Range<Position>,
    pub cnts: Stranded<Option<&'a [NucCounts]>>,
    pub coverage: Stranded<u32>,
}

pub struct NucCounterResult<'a, Data> {
    pub contig: &'a str,
    pub mapped: Stranded<u32>,
    pub cnts: Vec<InnerNucCounts<'a, Data>>,
}
