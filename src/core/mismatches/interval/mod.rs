use std::ops::Range;

use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{mock, predicate::*};

pub use borrowed::RefIntervalMismatches;
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
}

pub trait SeparableIntervalMismatches: Sized {
    fn split(self, indices: &[usize]) -> Vec<Self>;
    fn extract(self, ranges: &[Range<usize>], into: &mut Vec<Self>);
}

pub trait IntermediateIntervalMismatches:
    IntervalMismatches + SeparableIntervalMismatches + IntermediateMismatches
{
}
impl<T: IntervalMismatches + SeparableIntervalMismatches + IntermediateMismatches> IntermediateIntervalMismatches
    for T
{
}

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
    }
    impl SeparableIntervalMismatches for IntervalMismatches {
        fn split(self, indices: &[usize]) -> Vec<MockIntervalMismatches>;
        fn extract(self, ranges: &[Range<usize>], into: &mut Vec<MockIntervalMismatches>);
    }
}
