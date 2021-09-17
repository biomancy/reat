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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ROI {
    pub interval: Interval,
    pub name: String,
}

impl PartialOrd<Self> for ROI {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ROI {
    fn cmp(&self, other: &Self) -> Ordering {
        let from_contig = self.interval.contig().cmp(other.interval.contig());
        if !from_contig.is_eq() {
            return from_contig;
        }

        let from_start = self.interval.range().start.cmp(&other.interval.range().start);
        if !from_start.is_eq() {
            return from_start;
        }

        let len = self.interval.range().end - self.interval.range().start;
        let otherlen = other.interval.range().end - other.interval.range().start;

        // ROI with larget size must go first!
        len.cmp(&otherlen)
    }
}

#[derive(Clone, PartialEq, Debug, Getters, Dissolve)]
pub struct ROIWorkload {
    interval: Interval,
    rois: Vec<ROI>,
}

impl ROIWorkload {
    pub fn new(interval: Interval, rois: Vec<ROI>) -> ROIWorkload {
        debug_assert!(rois.iter().map(|x| &x.interval).all(|x| {
            x.contig() == interval.contig()
                && x.range().start >= interval.range().start
                && x.range().end <= interval.range().end
        }));
        ROIWorkload { interval: interval, rois }
    }

    fn infer_buffer_bin(seed: &Interval, binsize: u32) -> Interval {
        let end = seed.range().end.max(seed.range().start + binsize as u64);
        Interval::new(seed.contig().to_string(), seed.range().start..end)
    }

    pub fn from_bed(bed: &Path, binsize: u32) -> Vec<ROIWorkload> {
        assert!(binsize > 0, "Binsize must be > 0");

        let extensions = bed.file_name().and_then(OsStr::to_str).expect(
            "Failed to infer extension for BED file with target regions. \
                          Only \".bed\" and \".bed.gz\" files are supported",
        );
        let extensions: Vec<&str> = extensions.split('.').collect();
        assert!(!extensions.is_empty());

        let file = File::open(bed).expect("Failed to open bed file");
        let file = io::BufReader::new(file);
        let last = *extensions.last().unwrap();

        if last == "gz" || last == "gzip" {
            assert_eq!(extensions[extensions.len() - 2], "bed");
            ROIWorkload::_from_bed(io::BufReader::new(GzDecoder::new(file)), binsize)
        } else {
            assert_eq!(last, "bed");
            ROIWorkload::_from_bed(file, binsize)
        }
    }

    fn _from_bed<T: BufRead>(mut reader: T, binsize: u32) -> Vec<ROIWorkload> {
        let mut rois = Vec::new();

        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to read BED file") != 0 {
            let line = buf.trim_end();
            if line.is_empty() {
                buf.clear();
                continue;
            }
            let split: Vec<&str> = line.split('\t').take(4).collect();
            assert!(split.len() >= 3);

            let start = split[1].parse().expect("Failed to reads string start");
            let end = split[2].parse().expect("Failed to reads string start");
            let interval = Interval::new(split[0].to_owned(), Range { start, end });

            let name = if split.len() >= 4 { split[3] } else { &"." };
            let name = name.to_string();

            rois.push(ROI { interval, name });
            buf.clear();
        }

        rois.sort();
        if rois.is_empty() {
            return Vec::new();
        }

        let mut workloads: Vec<ROIWorkload> = Vec::new();

        let mut bin = ROIWorkload::infer_buffer_bin(&rois[0].interval, binsize);
        let mut buffer = vec![rois[0].clone()];
        let mut maxend = buffer[0].interval.range().end;
        for roi in rois.into_iter().skip(1) {
            // outside of the bin
            if roi.interval.contig() != bin.contig() || roi.interval.range().end > bin.range().end {
                debug_assert!(!buffer.is_empty());
                if maxend < bin.range().end {
                    bin.range_mut().end = maxend
                }

                workloads.push(ROIWorkload::new(bin, buffer));
                bin = ROIWorkload::infer_buffer_bin(&roi.interval, binsize);
                buffer = vec![roi];
                maxend = buffer[0].interval.range().end;
            } else {
                maxend = maxend.max(roi.interval.range().end);
                buffer.push(roi)
            }
        }

        // Save results for the last bin
        if maxend < bin.range().end {
            bin.range_mut().end = maxend
        }
        workloads.push(ROIWorkload::new(bin, buffer));

        workloads
    }

    fn _bin_chromosome(chromosome: &str, len: u64, binsize: u64) -> Vec<Interval> {
        let total_bins = (len + binsize - 1) / binsize;
        let mut bins: Vec<Interval> = Vec::with_capacity(total_bins as usize);

        for bin in 1..(total_bins + 1) {
            let start = (bin - 1) * binsize;
            let end = min(bin * binsize, len);
            debug_assert!(end - start <= binsize);

            bins.push(Interval::new(chromosome.to_owned(), start..end))
        }
        bins
    }

    pub fn from_hts(hts: &[impl AsRef<Path>], binsize: u64) -> Vec<ROIWorkload> {
        assert!(binsize > 0, "Binsize must be > 0");

        let hts = hts
            .iter()
            .map(|x| {
                IndexedReader::from_path(x).unwrap_or_else(|_| panic!("Failed to open file {}", x.as_ref().display()))
            })
            .collect_vec();

        // Check that headers are identical
        let all_equal = hts
            .iter()
            .map(|x| x.header())
            .map(|h| {
                (0..h.target_count())
                    .map(|tid| (String::from_utf8_lossy(h.tid2name(tid)), h.target_len(tid)))
                    .sorted()
                    .collect_vec()
            })
            .all_equal();
        assert!(all_equal, "BAM files must be mapped against identical reference assemblies");

        let header = hts[0].header();
        let mut workloads = Vec::with_capacity((4f32 * 10f32.powi(10) / binsize as f32) as usize);
        for tid in 0..header.target_count() {
            let tname = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            let tlen = header.target_len(tid).unwrap();
            workloads.extend(
                ROIWorkload::_bin_chromosome(&tname, tlen, binsize)
                    .into_iter()
                    .map(|e| ROIWorkload { interval: e.clone(), rois: vec![ROI { interval: e, name: "".into() }] }),
            )
        }
        workloads
    }

    pub fn len(&self) -> u64 {
        self.interval.range().end - self.interval.range().start
    }
}

impl Default for ROIWorkload {
    fn default() -> Self {
        Self { interval: Interval::new("".into(), 0..0), rois: vec![] }
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use bio_types::genome::Position;

    use super::*;

    fn interval(chr: &str, range: Range<Position>) -> Interval {
        Interval::new(chr.to_string(), range)
    }

    fn roi(chr: &str, range: Range<Position>, name: &str) -> ROI {
        ROI { interval: interval(chr, range), name: name.to_string() }
    }

    fn validate(bed: &str, expected: Vec<ROIWorkload>, binsizes: &[u32]) {
        for binsize in binsizes {
            let result = ROIWorkload::_from_bed(BufReader::new(bed.as_bytes()), *binsize);
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn non_overlapping() {
        let bed = "\
        chr1\t10\t20\tReg1\n\
        chr1\t50\t60\tIII\n\
        chr1\t30\t40\t2\n\
        chr1\t70\t80\t\n";
        let rois = vec![
            roi("chr1", 10..20, "Reg1"),
            roi("chr1", 30..40, "2"),
            roi("chr1", 50..60, "III"),
            roi("chr1", 70..80, "."),
        ];

        let expected = vec![
            ROIWorkload::new(interval("chr1", 10..20), vec![rois[0].clone()]),
            ROIWorkload::new(interval("chr1", 30..40), vec![rois[1].clone()]),
            ROIWorkload::new(interval("chr1", 50..60), vec![rois[2].clone()]),
            ROIWorkload::new(interval("chr1", 70..80), vec![rois[3].clone()]),
        ];
        validate(bed, expected, &[1, 15, 29]);

        let expected = vec![
            ROIWorkload::new(interval("chr1", 10..40), rois[0..2].to_vec()),
            ROIWorkload::new(interval("chr1", 50..80), rois[2..].to_vec()),
        ];
        validate(bed, expected, &[30, 40, 49]);

        let expected = vec![ROIWorkload::new(interval("chr1", 10..80), rois)];
        validate(bed, expected, &[70]);
    }

    #[test]
    fn overlapping() {
        let bed = "\
        1\t0\t3\t-\n\
        1\t3\t7\t+\n\
        1\t2\t5\t-\n\
        1\t4\t8\t+\n\
        1\t4\t8\t+\n";
        let rois = vec![
            roi("1", 0..3, "-"),
            roi("1", 2..5, "-"),
            roi("1", 3..7, "+"),
            roi("1", 4..8, "+"),
            roi("1", 4..8, "+"),
        ];

        let expected = vec![
            ROIWorkload::new(interval("1", 0..3), vec![rois[0].clone()]),
            ROIWorkload::new(interval("1", 2..5), vec![rois[1].clone()]),
            ROIWorkload::new(interval("1", 3..7), vec![rois[2].clone()]),
            ROIWorkload::new(interval("1", 4..8), vec![rois[3].clone(), rois[4].clone()]),
        ];
        validate(bed, expected, &[1, 2, 4]);

        let expected = vec![
            ROIWorkload::new(interval("1", 0..5), rois[0..2].to_vec()),
            ROIWorkload::new(interval("1", 3..8), rois[2..].to_vec()),
        ];
        validate(bed, expected, &[5, 6]);

        let expected = vec![ROIWorkload::new(interval("1", 0..8), rois)];
        validate(bed, expected, &[8, 100]);
    }

    #[test]
    fn none() {
        let bed = "\n\n\n";

        let expected: Vec<ROIWorkload> = vec![];
        validate(bed, expected, &[1, 2, 4]);
    }

    #[test]
    fn single() {
        let bed = "\nMT\t100\t200\t\n\n\n";
        let expected = vec![ROIWorkload::new(interval("MT", 100..200), vec![roi("MT", 100..200, ".")])];

        validate(bed, expected, &[1, 100, 500]);
    }

    #[test]
    fn complex() {
        let bed = "\
        1\t1\t3\t1-1\n\
        1\t2\t3\t1-2\n\
        1\t100\t110\t1-3\n\
        \n\
        2\t110\t120\t\n\
        2\t0\t200\t\n\
        2\t30\t50\t\n\
        \n\
        3\t10\t20\t2\n\
        3\t20\t30\t3\n\
        3\t0\t10\t1\n\
        3\t30\t40\t4\n\n\n";
        let rois = vec![
            roi("1", 1..3, "1-1"),
            roi("1", 2..3, "1-2"),
            roi("1", 100..110, "1-3"),
            roi("2", 0..200, "."),
            roi("2", 30..50, "."),
            roi("2", 110..120, "."),
            roi("3", 0..10, "1"),
            roi("3", 10..20, "2"),
            roi("3", 20..30, "3"),
            roi("3", 30..40, "4"),
        ];

        let expected = vec![
            ROIWorkload::new(interval("1", 1..110), rois[..3].to_vec()),
            ROIWorkload::new(interval("2", 0..200), rois[3..6].to_vec()),
            ROIWorkload::new(interval("3", 0..40), rois[6..].to_vec()),
        ];
        validate(bed, expected, &[200, 300]);

        let expected = vec![
            ROIWorkload::new(interval("1", 1..3), rois[..2].to_vec()),
            ROIWorkload::new(interval("1", 100..110), rois[2..3].to_vec()),
            ROIWorkload::new(interval("2", 0..200), rois[3..6].to_vec()),
            ROIWorkload::new(interval("3", 0..40), rois[6..].to_vec()),
        ];
        validate(bed, expected, &[50]);
    }

    #[test]
    fn bin_chromosome() {
        let workload = ROIWorkload::_bin_chromosome("chr1", 284, 100);
        let expected = vec![interval("chr1", 0..100), interval("chr1", 100..200), interval("chr1", 200..284)];

        assert_eq!(workload, expected);

        let workload = ROIWorkload::_bin_chromosome("2", 10, 1000);
        assert_eq!(workload.len(), 1);
        assert_eq!(workload[0], Interval::new("2".into(), 0..10));
    }
}
