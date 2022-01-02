use bio_types::genome::{AbstractInterval, Interval, Locus, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::Itertools;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::interval::BorrowedIntervalMismatches;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{NucCounter, NucCountingResult};
use crate::core::rpileup::ReadsCollider;

use super::base::BaseNucCounter;

pub struct IntervalNucCounts<'a> {
    block: &'a Interval,
    ncounts: &'a [NucCounts],
    strand: Strand,
}

impl<'a> NucCountingResult<'a, BorrowedIntervalMismatches<'a>> for IntervalNucCounts<'a> {
    fn interval(&self) -> &Interval {
        &self.block
    }

    fn ncounts(&self) -> &[NucCounts] {
        &self.ncounts
    }

    fn mismatches(self, reference: &'a [Nucleotide]) -> Vec<BorrowedIntervalMismatches<'a>> {
        let object = BorrowedIntervalMismatches::new(self.block.clone(), self.strand, reference, self.ncounts);
        vec![object]
    }
}

pub struct IntervalNucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for IntervalNucCounter<R, Filter>
where
    Self: 'a,
{
    type ColliderResult = IntervalNucCounts<'a>;
    type Workload = Interval;

    fn reset(&mut self, info: Self::Workload) {
        self.base.reset(info);
    }

    fn collide(&mut self, read: &R) {
        self.base.count(read);
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Vec<Self::ColliderResult> {
        vec![IntervalNucCounts { block: self.base.interval(), ncounts: self.base.counted(), strand: Strand::Unknown }]
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R, BorrowedIntervalMismatches<'a>>
    for IntervalNucCounter<R, Filter>
{
}
