use bio_types::genome::Interval;
use bio_types::strand::{ReqStrand, Strand};
use derive_getters::{Dissolve, Getters};

pub use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;

pub mod cnt;
pub mod filters;

pub trait NucCounter<'a, R: AlignedRead, Mismatches>: ReadsCollider<'a, R>
where
    Self::ColliderResult: NucCountingResult<'a, Mismatches>,
{
}

pub trait NucCountingResult<'a, Mismatches> {
    fn interval(&self) -> &Interval;
    fn ncounts(&self) -> &[NucCounts];
    fn mismatches(self, reference: &'a [Nucleotide]) -> Vec<Mismatches>;
}
