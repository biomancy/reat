use bio_types::genome::{AbstractInterval, Interval, Locus, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::Itertools;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::RefROIMismatches;
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{NucCounter, NucCountingResult};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::{ROIWorkload, ROI};

use super::base::BaseNucCounter;

pub struct ROINucCounts<'a> {
    coverage: u32,
    interval: &'a Interval,
    name: &'a String,
    strand: Strand,
    ncounts: &'a [NucCounts],
}

impl<'a> NucCountingResult<'a, RefROIMismatches<'a>> for ROINucCounts<'a> {
    fn interval(&self) -> &Interval {
        self.interval()
    }

    fn ncounts(&self) -> &[NucCounts] {
        &self.ncounts
    }

    fn mismatches(self, reference: &'a [Nucleotide]) -> Vec<RefROIMismatches<'a>> {
        vec![RefROIMismatches::new(self.interval, self.strand, self.name, self.coverage, reference, self.ncounts)]
    }
}

pub struct ROINucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    rois: Vec<ROI>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>> ROINucCounter<R, Filter> {
    pub fn new(maxbuf: usize, filter: Filter, trim5: u16, trim3: u16) -> Self {
        Self { base: BaseNucCounter::new(maxbuf, filter, trim5, trim3), rois: vec![] }
    }
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for ROINucCounter<R, Filter> {
    type ColliderResult = ROINucCounts<'a>;
    type Workload = ROIWorkload;

    fn reset(&mut self, info: Self::Workload) {
        let (interval, rois) = info.dissolve();
        self.base.reset(interval);
        self.rois = rois;
    }

    fn collide(&mut self, read: &R) {
        self.base.count(read);
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Vec<Self::ColliderResult> {
        let ncounts = self.base.counted();
        let interval = self.base.interval();
        let bufstart = interval.range().start as usize;

        let mut results = Vec::with_capacity(self.rois.len());
        for roi in &self.rois {
            debug_assert_eq!(roi.interval.contig(), interval.contig());
            let (start, end) = (roi.interval.range().start as usize, roi.interval.range().end as usize);
            let roicnts = &ncounts[start - bufstart..end - bufstart];
            // TODO: Add coverage values to this stuff
            todo!();
            results.push(ROINucCounts {
                coverage: 0,
                interval: &roi.interval,
                name: &roi.name,
                strand: Strand::Unknown,
                ncounts: roicnts,
            });
        }
        results
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R, RefROIMismatches<'a>>
    for ROINucCounter<R, Filter>
{
}
