use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Locus, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::Itertools;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::interval::RefIntervalMismatches;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{CountingResults, GroupedNucCounts, NucCounter, ToMismatches};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::SiteWorkload;

use super::base::BaseNucCounter;

pub struct IntervalNucCounts<'a> {
    block: Interval,
    ncounts: &'a [NucCounts],
    strand: Strand,
}

impl<'a> AbstractInterval for IntervalNucCounts<'a> {
    fn contig(&self) -> &str {
        self.block.contig()
    }

    fn range(&self) -> Range<Position> {
        self.block.range()
    }
}

impl<'a> CountingResults<'a> for IntervalNucCounts<'a> {
    fn ncounts(&self) -> &[NucCounts] {
        &self.ncounts
    }
}

impl<'a> ToMismatches<'a, RefIntervalMismatches<'a>> for IntervalNucCounts<'a> {
    fn mismatches(self, reference: &'a [Nucleotide]) -> Vec<RefIntervalMismatches<'a>> {
        vec![RefIntervalMismatches::new(self.block, self.strand, reference, self.ncounts)]
    }
}

pub struct IntervalNucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    ranges: Vec<Range<u64>>,
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for IntervalNucCounter<R, Filter>
where
    Self: 'a,
{
    type ColliderResult = GroupedNucCounts<'a, IntervalNucCounts<'a>>;
    type Workload = SiteWorkload;

    fn reset(&mut self, info: Self::Workload) {
        let (interval, ranges) = info.dissolve();
        self.base.reset(interval);
        self.ranges = ranges;
    }

    fn collide(&mut self, read: &R) {
        self.base.count(read);
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Self::ColliderResult {
        let (contig, start) = (self.base.interval().contig(), self.base.interval().range().start);
        let items = self
            .ranges
            .iter()
            .map(|x| {
                let indx = (x.start - start) as usize..(x.end - start) as usize;
                IntervalNucCounts {
                    block: Interval::new(contig.into(), x.clone()),
                    ncounts: &self.base.counted()[indx],
                    strand: Strand::Unknown,
                }
            })
            .collect();
        GroupedNucCounts::new(self.base.interval().clone(), items)
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R> for IntervalNucCounter<R, Filter> {
    type NucCounts = IntervalNucCounts<'a>;
}
