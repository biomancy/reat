



use std::io::Write;



use crate::core::strandutil::Stranded;

pub mod prefilters;
pub mod roi;
pub mod site;

pub type StrandingCounts = Stranded<usize>;

pub trait MismatchesVec: Sized {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;

    fn to_csv<F: Write>(&self, writer: &mut csv::Writer<F>) -> csv::Result<()>;
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
