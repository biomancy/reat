use std::cmp::Ordering;
use std::cmp::Ordering::Greater;
use std::cmp::{max, min};
use std::collections::{BTreeMap, HashMap};
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

#[derive(Eq, PartialEq, Debug)]
pub struct Bin<T> {
    pub bin: Interval,
    pub items: Vec<T>,
}

fn cmp(this: &impl AbstractInterval, other: &impl AbstractInterval) -> Ordering {
    let from_contig = this.contig().cmp(other.contig());
    if !from_contig.is_eq() {
        return from_contig;
    }

    let from_start = this.range().start.cmp(&other.range().start);
    if !from_start.is_eq() {
        return from_start;
    }

    let len = this.range().end - this.range().start;
    let otherlen = other.range().end - other.range().start;

    // Intervals with larget size must go first!
    len.cmp(&otherlen)
}

fn infer_bin<T: AbstractInterval>(seed: &T, binsize: u64) -> Interval {
    let end = seed.range().end.max(seed.range().start + binsize);
    Interval::new(seed.contig().to_string(), seed.range().start..end)
}

pub fn bin<T: AbstractInterval>(mut workloads: Vec<T>, maxbinsize: u64) -> Vec<Bin<T>> {
    if workloads.is_empty() {
        return vec![];
    }
    workloads.sort_by(cmp);

    let mut result = Vec::with_capacity(workloads.len());

    let mut iter = workloads.into_iter();
    let first = iter.next().unwrap();
    let mut bin = infer_bin(&first, maxbinsize);
    let mut buffer = vec![first];
    let mut maxend = buffer[0].range().end;

    for work in iter {
        // outside of the bin
        if work.contig() != bin.contig() || work.range().end > bin.range().end {
            debug_assert!(!buffer.is_empty());
            if maxend < bin.range().end {
                bin.range_mut().end = maxend
            }

            result.push(Bin { bin, items: buffer });
            bin = infer_bin(&work, maxbinsize);
            buffer = vec![work];
            maxend = buffer[0].range().end;
        } else {
            maxend = maxend.max(work.range().end);
            buffer.push(work)
        }
    }

    // Save results for the last bin
    if maxend < bin.range().end {
        bin.range_mut().end = maxend
    }
    result.push(Bin { bin, items: buffer });

    result
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use bio_types::genome::Position;

    use super::*;

    fn interval(chr: &str, range: Range<Position>) -> Interval {
        Interval::new(chr.to_string(), range)
    }

    fn validate(input: Vec<Interval>, expected: Vec<Bin<Interval>>, binsizes: &[u64]) {
        for binsize in binsizes {
            let result = bin(input.clone(), *binsize);
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn non_overlapping() {
        let inter = vec![
            interval("chr1", 10..20),
            interval("chr1", 30..40),
            interval("chr1", 50..60),
            interval("chr1", 70..80),
        ];

        let expected = vec![
            Bin { bin: interval("chr1", 10..20), items: vec![inter[0].clone()] },
            Bin { bin: interval("chr1", 30..40), items: vec![inter[1].clone()] },
            Bin { bin: interval("chr1", 50..60), items: vec![inter[2].clone()] },
            Bin { bin: interval("chr1", 70..80), items: vec![inter[3].clone()] },
        ];
        validate(inter.clone(), expected, &[1, 15, 29]);

        let expected = vec![
            Bin { bin: interval("chr1", 10..40), items: inter[0..2].to_vec() },
            Bin { bin: interval("chr1", 50..80), items: inter[2..].to_vec() },
        ];
        validate(inter.clone(), expected, &[30, 40, 49]);

        let expected = vec![Bin { bin: interval("chr1", 10..80), items: inter.clone() }];
        validate(inter, expected, &[70]);
    }

    #[test]
    fn overlapping() {
        let inter = vec![
            interval("1", 0..3),
            interval("1", 2..5),
            interval("1", 3..7),
            interval("1", 4..8),
            interval("1", 4..8),
        ];

        let expected = vec![
            Bin { bin: interval("1", 0..3), items: vec![inter[0].clone()] },
            Bin { bin: interval("1", 2..5), items: vec![inter[1].clone()] },
            Bin { bin: interval("1", 3..7), items: vec![inter[2].clone()] },
            Bin { bin: interval("1", 4..8), items: vec![inter[3].clone(), inter[4].clone()] },
        ];
        validate(inter.clone(), expected, &[1, 2, 4]);

        let expected = vec![
            Bin { bin: interval("1", 0..5), items: inter[0..2].to_vec() },
            Bin { bin: interval("1", 3..8), items: inter[2..].to_vec() },
        ];
        validate(inter.clone(), expected, &[5, 6]);

        let expected = vec![Bin { bin: interval("1", 0..8), items: inter.clone() }];
        validate(inter, expected, &[8, 100]);
    }

    #[test]
    fn none() {
        validate(vec![], vec![], &[1, 2, 3, 4]);
    }

    #[test]
    fn single() {
        let inter = interval("MT", 100..200);
        let expected = vec![Bin { bin: interval("MT", 100..200), items: vec![inter.clone()] }];

        validate(vec![inter], expected, &[1, 100, 500]);
    }

    #[test]
    fn complex() {
        let inter = vec![
            interval("1", 1..3),
            interval("1", 2..3),
            interval("1", 100..110),
            interval("2", 0..200),
            interval("2", 30..50),
            interval("2", 110..120),
            interval("3", 0..10),
            interval("3", 10..20),
            interval("3", 20..30),
            interval("3", 30..40),
        ];

        let expected = vec![
            Bin { bin: interval("1", 1..110), items: inter[..3].to_vec() },
            Bin { bin: interval("2", 0..200), items: inter[3..6].to_vec() },
            Bin { bin: interval("3", 0..40), items: inter[6..].to_vec() },
        ];
        validate(inter.clone(), expected, &[200, 300]);

        let expected = vec![
            Bin { bin: interval("1", 1..3), items: inter[..2].to_vec() },
            Bin { bin: interval("1", 100..110), items: inter[2..3].to_vec() },
            Bin { bin: interval("2", 0..200), items: inter[3..6].to_vec() },
            Bin { bin: interval("3", 0..40), items: inter[6..].to_vec() },
        ];
        validate(inter, expected, &[50]);
    }
}
