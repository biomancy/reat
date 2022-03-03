use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};

use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, AggregatedNucCountsItem, NucCounter};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::SiteWorkload;

use super::base::BaseNucCounter;

pub struct IntervalNucCounts<'a> {
    contig: &'a str,
    range: Range<Position>,
    items: Vec<AggregatedNucCountsItem<'a, ()>>,
}

impl<'a> AggregatedNucCounts<'a> for IntervalNucCounts<'a> {
    type ItemInfo = ();

    fn items(&self) -> &[AggregatedNucCountsItem<'a, Self::ItemInfo>] {
        &self.items
    }

    fn items_mut(&mut self) -> &mut [AggregatedNucCountsItem<'a, Self::ItemInfo>] {
        &mut self.items
    }

    fn consume(self) -> (&'a str, Vec<AggregatedNucCountsItem<'a, Self::ItemInfo>>) {
        (self.contig, self.items)
    }
}

impl<'a> AbstractInterval for IntervalNucCounts<'a> {
    fn contig(&self) -> &str {
        self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

#[derive(Clone)]
pub struct IntervalNucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    ranges: Vec<Range<u64>>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>> IntervalNucCounter<R, Filter> {
    pub fn new(base: BaseNucCounter<R, Filter>) -> Self {
        Self { base, ranges: vec![] }
    }
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for IntervalNucCounter<R, Filter> {
    type ColliderResult = IntervalNucCounts<'a>;
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
        let (contig, range) = (self.base.interval().contig(), self.base.interval().range());
        let start = range.start;

        let items = self
            .ranges
            .iter()
            .map(|x| {
                let indx = (x.start - start) as usize..(x.end - start) as usize;
                AggregatedNucCountsItem {
                    info: (),
                    range: x.clone(),
                    forward: None,
                    reverse: None,
                    unstranded: Some(&self.base.counted()[indx]),
                }
            })
            .collect();
        Self::ColliderResult { contig, range, items }
    }
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> NucCounter<'a, R> for IntervalNucCounter<R, Filter> {}
