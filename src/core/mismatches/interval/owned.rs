use core::borrow::Borrow;

use bio_types::genome::Interval;
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::interval::borrowed::BorrowedIntervalMismatches;
use crate::core::mismatches::interval::IntervalMismatches;

#[derive(Constructor)]
pub struct OwnedIntervalMismatches {
    interval: Interval,
    strand: Strand,
    refnuc: Vec<Nucleotide>,
    ncounts: Vec<NucCounts>,
}

impl IntervalMismatches for OwnedIntervalMismatches {
    fn interval(&self) -> &Interval {
        &self.interval
    }

    fn strand(&self) -> &Strand {
        &self.strand
    }

    fn refnuc(&self) -> &[Nucleotide] {
        &self.refnuc
    }

    fn ncounts(&self) -> &[NucCounts] {
        &self.ncounts
    }

    fn split(self, indices: &[usize]) -> Vec<Self>
    where
        Self: Sized,
    {
        todo!()
    }
}

impl From<BorrowedIntervalMismatches<'_>> for OwnedIntervalMismatches {
    fn from(x: BorrowedIntervalMismatches<'_>) -> Self {
        OwnedIntervalMismatches::new(
            x.interval().to_owned(),
            x.strand().clone(),
            x.refnuc().to_owned(),
            x.ncounts().to_owned(),
        )
    }
}
