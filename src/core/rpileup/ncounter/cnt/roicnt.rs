use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use itertools::{zip, Itertools};
use serde::de::Unexpected::Str;

use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounter::filters::ReadsFilter;
use crate::core::rpileup::ncounter::{InnerNucCounts, NucCounterResult};
use crate::core::rpileup::ReadsCollider;
use crate::core::strandutil::Stranded;
use crate::core::workload::{ROIWorkload, ROI};

use super::base::BaseNucCounter;

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
    type ColliderResult = NucCounterResult<'a, &'a ROI>;
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
        let contig = self.base.interval().contig();
        let instart = self.base.interval().range().start as usize;

        let mut cnts = Vec::with_capacity(self.rois.len());
        for (coverage, roi) in zip(&self.coverage, &self.rois) {
            debug_assert_eq!(roi.contig(), contig);
            let (start, end) = (roi.range().start as usize, roi.range().end as usize);

            let roicnts = &self.base.counted()[start - instart..end - instart];
            cnts.push(InnerNucCounts {
                data: roi,
                range: roi.range().clone(),
                cnts: Stranded::unknown(Some(roicnts)),
                coverage: Stranded::unknown(*coverage),
            });
        }
        NucCounterResult { contig, cnts, mapped: Stranded::unknown(self.base.mapped()) }
    }
}
