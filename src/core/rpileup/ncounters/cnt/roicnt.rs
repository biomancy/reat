use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Position};
use itertools::{zip, Itertools};

use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, AggregatedNucCountsItem, NucCounter};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::{ROIWorkload, ROI};

use super::base::BaseNucCounter;
use crate::core::strandutil::StrandedData;

#[derive(PartialEq)]
pub struct ROINucCountsInfo<'a> {
    pub coverage: u32,
    pub roi: &'a ROI,
}

pub struct ROINucCounts<'a> {
    contig: &'a str,
    range: Range<Position>,
    items: Vec<AggregatedNucCountsItem<'a, ROINucCountsInfo<'a>>>,
}

impl<'a> AbstractInterval for ROINucCounts<'a> {
    fn contig(&self) -> &str {
        self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

impl<'a> AggregatedNucCounts<'a> for ROINucCounts<'a> {
    type ItemInfo = ROINucCountsInfo<'a>;

    fn items(&self) -> &[AggregatedNucCountsItem<'a, Self::ItemInfo>] {
        &self.items
    }

    fn items_mut(&mut self) -> &mut [AggregatedNucCountsItem<'a, Self::ItemInfo>] {
        &mut self.items
    }

    fn consume(self) -> (&'a str, Vec<AggregatedNucCountsItem<'a, Self::ItemInfo>>) {
        (&self.contig, self.items)
    }
}

#[derive(Clone)]
pub struct ROINucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    rois: Vec<ROI>,
    coverage: Vec<u32>,
    index: IntervalTree<u32, usize>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>> ROINucCounter<R, Filter> {
    pub fn new(base: BaseNucCounter<R, Filter>) -> Self {
        Self { base, rois: vec![], coverage: vec![], index: Default::default() }
    }
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for ROINucCounter<R, Filter> {
    type ColliderResult = ROINucCounts<'a>;
    type Workload = ROIWorkload;

    fn reset(&mut self, info: Self::Workload) {
        let (bin, rois) = info.dissolve();
        self.base.reset(bin);
        self.rois = rois;

        // Coverage for each roi
        self.coverage.clear();
        self.coverage.resize(self.rois.len(), 0);

        // Index to accurately count ROIs coverage
        self.index = Default::default();
        let binstart = self.base.interval().range().start;
        for (ind, roi) in self.rois.iter().enumerate() {
            let (start, end) = (roi.range().start - binstart, roi.range().end - binstart);
            self.index.insert(start as u32..end as u32, ind)
        }
    }

    fn collide(&mut self, read: &R) {
        let covered_rois =
            self.base.count(read).iter().map(|x| self.index.find(x)).flatten().map(|x| x.data()).unique();
        for ind in covered_rois {
            self.coverage[*ind] += 1;
        }
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Self::ColliderResult {
        let (contig, range) = (self.base.interval().contig(), self.base.interval().range());
        let instart = range.start as usize;

        let mut items = Vec::with_capacity(self.rois.len());
        for (coverage, roi) in zip(&self.coverage, &self.rois) {
            debug_assert_eq!(roi.contig(), contig);
            let info = ROINucCountsInfo { coverage: *coverage, roi };
            let (start, end) = (roi.range().start as usize, roi.range().end as usize);
            let roicnts = &self.base.counted()[start - instart..end - instart];

            items.push(AggregatedNucCountsItem {
                info,
                range: roi.range(),
                mapped: StrandedData { forward: 0, reverse: 0, unknown: self.base.mapped() },
                counts: StrandedData { forward: None, reverse: None, unknown: Some(roicnts) },
            })
        }
        ROINucCounts { contig, range, items }
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R> for ROINucCounter<R, Filter> {}
