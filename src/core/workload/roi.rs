use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{Same, Strand};
use derive_getters::{Dissolve, Getters};

use crate::core::io::bed::BedRecord;

use super::utils;

#[derive(Clone, Debug, Getters, Dissolve)]
pub struct ROI {
    interval: Interval,
    name: String,
    strand: Strand,
    subintervals: Vec<Range<u64>>,
}

impl PartialEq for ROI {
    fn eq(&self, other: &Self) -> bool {
        self.interval == other.interval
            && self.name == other.name
            && self.strand.same(&other.strand)
            && self.subintervals == other.subintervals
    }
}

impl AbstractInterval for ROI {
    fn contig(&self) -> &str {
        self.interval.contig()
    }

    fn range(&self) -> Range<Position> {
        self.subintervals.first().unwrap().start..self.subintervals.last().unwrap().end
    }
}

impl ROI {
    pub fn new(roi: Interval, name: String, strand: Strand, subintervals: Vec<Range<u64>>) -> Self {
        debug_assert!(!subintervals.is_empty());
        debug_assert!(subintervals.iter().all(|x| x.start >= roi.range().start && x.end <= roi.range().end));
        ROI { interval: roi, name, strand, subintervals }
    }

    pub fn nucmasked(&self) -> u32 {
        let mut nucin = 0;
        for piece in &self.subintervals {
            nucin += piece.end - piece.start;
        }
        let total = self.interval.range().end - self.interval.range().start;

        (total - nucin) as u32
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
        debug_assert!(rois.iter().all(|x| bin.contig() == x.interval.contig()
            && bin.range().contains(&x.range().start)
            && bin.range().contains(&x.range().end)));
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
                .map(|x| ROI::new(x.inner.interval, x.inner.name, x.inner.strand, x.retained))
                .collect()
        } else {
            rois.into_iter()
                .map(|x| {
                    let range = x.interval.range();
                    ROI::new(x.interval, x.name, x.strand, vec![range])
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
