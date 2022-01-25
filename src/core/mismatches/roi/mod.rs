use std::borrow::Borrow;
use std::path::Path;

use bio_types::genome::Interval;
use bio_types::strand::{Same, Strand};
#[cfg(test)]
use mockall::{mock, predicate::*};

pub use borrowed::RefROIMismatches;
pub use msummary::MismatchesSummary;
pub use owned::OwnedROIMismatches;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::mismatches::IntermediateMismatches;
use crate::core::workload::ROI;

pub trait ROIMismatches {
    fn roi(&self) -> &ROI;
    fn strand(&self) -> &Strand;
    fn masked(&self) -> u32;
    fn coverage(&self) -> u32;
    fn sequence(&self) -> &NucCounts;
    fn mismatches(&self) -> &MismatchesSummary;
}

impl PartialEq for dyn ROIMismatches {
    fn eq(&self, other: &Self) -> bool {
        self.roi() == other.roi()
            && self.strand().same(other.strand())
            && self.masked() == other.masked()
            && self.coverage() == other.coverage()
            && self.sequence() == other.sequence()
            && self.mismatches() == other.mismatches()
    }
}

mod borrowed;
mod msummary;
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
        fn roi(&self) -> &ROI;
        fn strand(&self) -> &Strand;
        fn masked(&self) -> u32;
        fn coverage(&self) -> u32;
        fn sequence(&self) -> &NucCounts;
        fn mismatches(&self) -> &MismatchesSummary;
    }
}
