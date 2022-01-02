use std::borrow::Borrow;
use std::path::Path;

use bio_types::genome::Interval;
use bio_types::strand::{Same, Strand};
#[cfg(test)]
use mockall::{mock, predicate::*};

pub use borrowed::BorrowedROIMismatches;
pub use mismatches_summary::MismatchesSummary;
pub use owned::OwnedROIMismatches;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::mismatches::IntermediateMismatches;

pub trait ROIMismatches {
    fn interval(&self) -> &Interval;
    fn strand(&self) -> &Strand;
    fn name(&self) -> &String;
    fn coverage(&self) -> u32;
    fn sequence(&self) -> &NucCounts;
    fn mismatches(&self) -> &MismatchesSummary;
}

impl PartialEq for dyn ROIMismatches {
    fn eq(&self, other: &Self) -> bool {
        self.interval() == other.interval()
            && self.strand().same(other.strand())
            && self.name() == other.name()
            && self.coverage() == other.coverage()
            && self.sequence() == other.sequence()
            && self.mismatches() == other.mismatches()
    }
}

mod borrowed;
mod mismatches_summary;
mod owned;

pub trait IntermediateROIMismatches: ROIMismatches + IntermediateMismatches {}
impl<T: ROIMismatches + IntermediateMismatches> IntermediateROIMismatches for T {}
#[cfg(test)]
mock! {
    pub ROIMismatches {}
    impl IntermediateMismatches for ROIMismatches {
        type Finalized = ();

        fn finalize(self) -> <Self as IntermediateMismatches>::Finalized;
        fn set_strand(&mut self, strand: Strand);
    }
    impl ROIMismatches for ROIMismatches {
        fn interval(&self) -> &Interval;
        fn strand(&self) -> &Strand;
        fn name(&self) -> &String;
        fn coverage(&self) -> u32;
        fn sequence(&self) -> &NucCounts;
        fn mismatches(&self) -> &MismatchesSummary;
    }
}
