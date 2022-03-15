use std::ops::Range;

use bio_types::genome::AbstractInterval;

use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounter::filters::ReadsFilter;
use crate::core::rpileup::ncounter::{InnerNucCounts, NucCounterResult};
use crate::core::rpileup::ReadsCollider;
use crate::core::strandutil::Stranded;
use crate::core::workload::SiteWorkload;

use super::base::BaseNucCounter;

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
    type ColliderResult = NucCounterResult<'a, ()>;
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
        let contig = self.base.interval().contig();
        let start = self.base.interval().range().start;

        let cnts = self
            .ranges
            .iter()
            .map(|range| {
                let indx = (range.start - start) as usize..(range.end - start) as usize;
                InnerNucCounts {
                    data: (),
                    range: range.clone(),
                    cnts: Stranded::unknown(Some(&self.base.counted()[indx])),
                    coverage: Stranded::unknown(self.base.mapped()),
                }
            })
            .collect();
        Self::ColliderResult { contig, cnts, mapped: Stranded::unknown(self.base.mapped()) }
    }
}
