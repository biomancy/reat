use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{Same, Strand};
use derive_getters::{Dissolve, Getters};
use soa_derive::StructOfArray;

use crate::core::io::bed::BedRecord;

use super::utils;

#[derive(Clone, Debug, Dissolve)]
pub struct ROI {
    contig: String,
    premasked: Range<Position>,
    subintervals: Vec<Range<Position>>,
    name: String,
    strand: Strand,
}

impl PartialEq for ROI {
    fn eq(&self, other: &Self) -> bool {
        self.contig == other.contig
            && self.premasked == other.premasked
            && self.strand.same(&other.strand)
            && self.name == other.name
            && self.subintervals == other.subintervals
    }
}

impl AbstractInterval for ROI {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.postmasked()
    }
}

impl ROI {
    pub fn new(
        contig: String,
        premasked: Range<Position>,
        subintervals: Vec<Range<Position>>,
        name: String,
        strand: Strand,
    ) -> Self {
        debug_assert!(!subintervals.is_empty());
        debug_assert!(subintervals.iter().all(|x| x.start >= premasked.start && x.end <= premasked.end));
        ROI { contig, premasked, subintervals, name, strand }
    }
    
    pub fn premasked(&self) -> Range<Position> {
        self.premasked.clone()
    }

    pub fn postmasked(&self) -> Range<Position> {
        self.subintervals.first().unwrap().start..self.subintervals.last().unwrap().end
    }

    pub fn subintervals(&self) -> &[Range<Position>] {
        &self.subintervals
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }
}

#[derive(Clone, PartialEq, Debug, Dissolve, Getters)]
pub struct ROIWorkload {
    bin: Interval,
    rois: Vec<ROI>,
}

impl AbstractInterval for ROIWorkload {
    fn contig(&self) -> &str {
        self.bin.contig()
    }

    fn range(&self) -> Range<Position> {
        self.bin.range()
    }
}

impl ROIWorkload {
    pub fn new(bin: Interval, rois: Vec<ROI>) -> Self {
        debug_assert!(rois.iter().all(|x| bin.contig() == x.contig()
            && bin.range().contains(&x.postmasked().start)
            && bin.range().contains(&x.postmasked().end)));
        ROIWorkload { bin, rois }
    }
}

#[allow(clippy::len_without_is_empty)]
impl ROIWorkload {
    pub fn from_bed(rois: Vec<BedRecord>, binsize: u64, exclude: Option<Vec<BedRecord>>) -> Vec<ROIWorkload> {
        assert!(binsize > 0, "Binsize must be > 0");

        // 1. Subtract from rois all the excluded regions and create ROI objects
        let rois = if let Some(exclude) = exclude {
            utils::subtract(rois, exclude)
                .into_iter()
                .map(|x| ROI::new(x.inner.contig().into(), x.inner.range(), x.retained, x.inner.name, x.inner.strand))
                .collect()
        } else {
            rois.into_iter()
                .map(|x| {
                    let subintervals = vec![x.interval.range()];
                    ROI::new(x.contig().into(), x.range(), subintervals, x.name, x.strand)
                })
                .collect()
        };

        // 2. Bin these guys and create workloads
        utils::bin(rois, binsize).into_iter().map(|x| ROIWorkload { bin: x.bin, rois: x.items }).collect()
    }

    #[inline]
    pub fn len(&self) -> usize {
        (self.bin.range().end - self.bin.range().start) as usize
    }
}
