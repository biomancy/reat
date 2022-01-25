use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{ReqStrand, Strand};
use derive_getters::{Dissolve, Getters};

pub use groupcnt::GroupedNucCounts;

pub use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ReadsCollider;

pub mod cnt;
pub mod filters;
mod groupcnt;

pub trait NucCounter<'a, R: AlignedRead>:
    ReadsCollider<'a, R, ColliderResult = GroupedNucCounts<'a, Self::NucCounts>>
{
    type NucCounts: CountingResults<'a>;
}

pub trait CountingResults<'a>: AbstractInterval {
    fn ncounts(&self) -> &[NucCounts];
}

pub trait ToMismatches<'a, T> {
    fn mismatches(self, reference: &'a [Nucleotide]) -> Vec<T>;
}
