use std::cmp::min;
use std::cmp::Ordering;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::ops::Range;
use std::path::Path;

use bio_types::genome::{AbstractInterval, Interval};
use derive_getters::{Dissolve, Getters};
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read};

use super::utils;

#[derive(Clone, PartialEq, Debug, Getters, Dissolve)]
pub struct SiteWorkload {
    interval: Interval,
    include: Vec<Range<u64>>,
}

impl SiteWorkload {
    pub fn from_intervals(
        mut intervals: Vec<Interval>,
        binsize: u64,
        exclude: Option<Vec<impl AbstractInterval>>,
    ) -> Vec<SiteWorkload> {
        assert!(binsize > 0, "Binsize must be > 0");
        // Subtract excluded if needed
        if let Some(bl) = exclude {
            intervals = utils::subtract(intervals, bl)
                .into_iter()
                .map(|x| {
                    let contig = x.inner.contig().to_owned();
                    x.retained.into_iter().map(move |piece| Interval::new(contig.clone(), piece))
                })
                .flatten()
                .collect();
        }

        // Bin and transform to the workload
        utils::bin(intervals, binsize)
            .into_iter()
            .map(|x| {
                let ranges = x.items.into_iter().map(|x| x.range()).collect();
                SiteWorkload { interval: x.bin, include: ranges }
            })
            .collect()
    }
}
