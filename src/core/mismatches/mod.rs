use bio_types::strand::Strand;
use std::io::Write;

use crate::core::strandutil::Stranded;

pub mod prefilters;
pub mod roi;
pub mod site;

pub type StrandingCounts = Stranded<usize>;

pub trait MismatchesVec: Sized {
    fn contig(&self) -> &str;
    fn trstrand(&self) -> Strand;

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;

    fn ugly_sort_and_to_csv<F: Write>(items: Vec<Self>, writer: &mut csv::Writer<F>) -> csv::Result<()>;
}

pub trait Builder<'a> {
    type Out: MismatchesVec;
    type SourceCounts;
    fn build(&mut self, nc: Self::SourceCounts) -> Batch<Self::Out>;
}

pub struct Batch<T: MismatchesVec> {
    pub contig: String,
    pub mapped: Stranded<u32>,
    // Must be retained & printed no matter what
    pub retained: Stranded<T>,
    // Other mismatches
    pub items: Stranded<T>,
}
