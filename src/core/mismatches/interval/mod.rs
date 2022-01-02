#[cfg(test)]
use mockall::{mock, predicate::*};

use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

pub use borrowed::BorrowedIntervalMismatches;
pub use owned::OwnedIntervalMismatches;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::mismatches::IntermediateMismatches;

mod borrowed;
mod owned;

pub trait IntervalMismatches {
    fn interval(&self) -> &Interval;
    fn strand(&self) -> &Strand;
    fn refnuc(&self) -> &[Nucleotide];
    fn ncounts(&self) -> &[NucCounts];

    fn split(self, indices: &[usize]) -> Vec<Self>
    where
        Self: Sized;
}

pub trait IntermediateIntervalMismatches: IntervalMismatches + IntermediateMismatches {}
impl<T: IntervalMismatches + IntermediateMismatches> IntermediateIntervalMismatches for T {}

#[cfg(test)]
mock! {
    pub IntervalMismatches {}
    impl IntermediateMismatches for IntervalMismatches {
        type Finalized = ();

        fn finalize(self) -> <Self as IntermediateMismatches>::Finalized;
        fn set_strand(&mut self, strand: Strand);
    }
    impl IntervalMismatches for IntervalMismatches {
        fn interval(&self) -> &Interval;
        fn strand(&self) -> &Strand;
        fn refnuc(&self) -> &[Nucleotide];
        fn ncounts(&self) -> &[NucCounts];

        fn split(self, indices: &[usize]) -> Vec<Self>
        where
            Self: Sized;
    }
}
