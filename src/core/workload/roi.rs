use std::cmp::min;
use std::cmp::Ordering;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::ops::Range;
use std::path::Path;

use bio_types::genome::{AbstractInterval, Interval, Position};
use derive_getters::{Dissolve, Getters};
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read};

use crate::core::io::bed;
use crate::core::io::bed::BedRecord;

use super::utils;

#[derive(Clone, Debug, Eq, PartialEq, Getters)]
pub struct ROI {
    include: Vec<Range<u64>>,
    original: Interval,
    name: String,
}

impl AbstractInterval for ROI {
    fn contig(&self) -> &str {
        self.original.contig()
    }

    fn range(&self) -> Range<Position> {
        self.include.first().unwrap().start..self.include.last().unwrap().end
    }
}

impl ROI {
    pub fn new(roi: Interval, name: String, include: Vec<Range<u64>>) -> Self {
        debug_assert!(!include.is_empty());
        debug_assert!(include.iter().all(|x| x.start >= roi.range().start && x.end <= roi.range().end));
        ROI { include, original: roi, name }
    }
}

#[derive(Clone, PartialEq, Debug, Dissolve, Getters)]
pub struct ROIWorkload {
    bin: Interval,
    rois: Vec<ROI>,
}

impl ROIWorkload {
    pub fn new(bin: Interval, rois: Vec<ROI>) -> Self {
        debug_assert!(rois.iter().all(|x| bin.contig() == x.original.contig()
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
        let rois = if let Some(bl) = exclude {
            utils::subtract(rois, bl)
                .into_iter()
                .map(|x| ROI::new(x.inner.interval, x.inner.name, x.retained))
                .collect()
        } else {
            rois.into_iter()
                .map(|x| {
                    let range = x.interval.range();
                    ROI::new(x.interval, x.name, vec![range])
                })
                .collect()
        };

        // 2. Bin these guys and create workloads
        utils::bin(rois, binsize).into_iter().map(|x| ROIWorkload { bin: x.bin, rois: x.items }).collect()
    }

    #[inline]
    pub fn len(&self) -> u64 {
        self.bin.range().end - self.bin.range().start
    }
}
