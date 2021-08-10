use std::cmp::min;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::ops::Range;
use std::path::{Path, PathBuf};

use bio_types::genome::{AbstractInterval, Interval};
use flate2::bufread::GzDecoder;
use rust_htslib::bam::{IndexedReader, Read};

#[derive(PartialEq, Debug)]
pub struct Workload {
    pub name: String,
    pub interval: Interval,
    pub bam: PathBuf,
}

impl Workload {
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

    pub fn from_binned_hts(hts: &Path, binsize: u64) -> Vec<Workload> {
        assert!(binsize > 0, "Binsize must be > 0");

        let bam = IndexedReader::from_path(hts).expect(&format!("Failed to open file {}", hts.display()));
        let header = bam.header();

        let mut workloads = Vec::with_capacity((4f32 * 10f32.powi(10) / binsize as f32) as usize);
        for tid in 0..header.target_count() {
            let tname = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            let tlen =
                header.target_len(tid).expect(&format!("Failed to parse {} length for file {}.", tname, hts.display()));

            workloads.extend(Workload::_bin_chromosome(&tname, tlen, binsize).into_iter().map(|e| Workload {
                name: format!("{}:{}-{}", e.contig(), e.range().start, e.range().start),
                interval: e,
                bam: hts.to_path_buf(),
            }))
        }
        workloads
    }

    fn _from_bed_intervals<T: BufRead>(mut reader: T, bam: &Path) -> Vec<Workload> {
        let mut result: Vec<Workload> = Vec::new();

        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to reads BED file") != 0 {
            let split: Vec<&str> = buf.trim_end().split('\t').take(4).collect();
            assert!(split.len() >= 3);

            if split[0].len() > 6 {
                buf.clear();
                continue;
            }

            let start = split[1].parse().expect("Failed to reads string start");
            let end = split[2].parse().expect("Failed to reads string start");
            let interval = Interval::new(split[0].to_owned(), Range { start, end });

            let name = if split.len() >= 4 { split[3] } else { &"." };

            result.push(Workload { name: name.to_string(), interval: interval, bam: bam.to_path_buf() });
            buf.clear();
        }
        result
    }

    pub fn from_bed_intervals<'a>(bed: &'a Path, bam: &Path) -> Vec<Workload> {
        let extensions = bed.file_name().and_then(OsStr::to_str).expect(
            "Failed to infer extension for BED file with target regions. \
                          Only \".bed\" and \".bed.gz\" files are supported",
        );
        let extensions: Vec<&str> = extensions.split('.').collect();
        assert!(extensions.len() >= 1);

        let file = File::open(bed).expect("Failed to open bed file");
        let file = io::BufReader::new(file);
        let last = *extensions.last().unwrap();

        if last == "gz" || last == "gzip" {
            assert_eq!(extensions[extensions.len() - 2], "bed");
            Workload::_from_bed_intervals(io::BufReader::new(GzDecoder::new(file)), bam)
        } else {
            assert_eq!(last, "bed");
            Workload::_from_bed_intervals(file, bam)
        }
    }

    pub fn len(&self) -> u64 {
        self.interval.range().end - self.interval.range().start
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;

    #[test]
    fn bin_chromosome() {
        let workload = Workload::_bin_chromosome("chr1", 284, 100);
        let expected = vec![
            Interval::new("chr1".into(), 0..100),
            Interval::new("chr1".into(), 100..200),
            Interval::new("chr1".into(), 200..284),
        ];

        assert_eq!(workload.len(), expected.len());
        for (work, exp) in workload.iter().zip(&expected) {
            assert_eq!(work, exp);
        }

        let workload = Workload::_bin_chromosome("2", 10, 1000);
        assert_eq!(workload.len(), 1);
        assert_eq!(workload[0], Interval::new("2".into(), 0..10));
    }

    #[test]
    fn from_bed_intervals() {
        let bed = "chr1\t100\t200\t1\n\
        chr2\t100\t200\tRegion 2\n\
        chr2\t150\t250\tRegion 3v.a\n\
        3\t1\t100\tVery long Region 4";
        let bam = PathBuf::from("test.bam");

        let workload = Workload::_from_bed_intervals(BufReader::new(bed.as_bytes()), &bam);

        let expected = vec![
            Workload { name: "1".into(), interval: Interval::new("chr1".into(), 100..200), bam: bam.clone() },
            Workload { name: "Region 2".into(), interval: Interval::new("chr2".into(), 100..200), bam: bam.clone() },
            Workload { name: "Region 3v.a".into(), interval: Interval::new("chr2".into(), 150..250), bam: bam.clone() },
            Workload {
                name: "Very long Region 4".into(),
                interval: Interval::new("3".into(), 1..100),
                bam: bam.clone(),
            },
        ];

        assert_eq!(workload.len(), expected.len());
        for (work, exp) in workload.iter().zip(&expected) {
            assert_eq!(work, exp);
        }
    }
}
