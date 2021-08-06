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

        let bam =
            IndexedReader::from_path(hts).expect(&format!("Failed to open file {}", hts.display()));
        let header = bam.header();

        let mut workloads = Vec::with_capacity((4f32 * 10f32.powi(10) / binsize as f32) as usize);
        for tid in 0..header.target_count() {
            let tname = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            let tlen = header.target_len(tid).expect(&format!(
                "Failed to parse {} length for file {}.",
                tname,
                hts.display()
            ));

            workloads.extend(
                Workload::_bin_chromosome(&tname, tlen, binsize)
                    .into_iter()
                    .map(|e| Workload {
                        name: format!("{}:{}-{}", e.contig(), e.range().start, e.range().start),
                        interval: e,
                        bam: hts.to_path_buf(),
                    }),
            )
        }
        workloads
    }

    fn _from_bed_interval<T: BufRead>(mut reader: T, bam: &Path) -> Vec<Workload> {
        let mut result: Vec<Workload> = Vec::new();

        let mut buf = String::new();
        while reader
            .read_line(&mut buf)
            .expect("Failed to records BED file")
            != 0
        {
            let split: Vec<&str> = buf.trim_end().split('\t').take(4).collect();
            assert!(split.len() >= 3);

            if split[0].len() > 6 {
                buf.clear();
                continue;
            }

            let start = split[1].parse().expect("Failed to records string start");
            let end = split[2].parse().expect("Failed to records string start");
            let interval = Interval::new(split[0].to_owned(), Range { start, end });

            let name = if split.len() >= 4 { split[3] } else { &"." };

            result.push(Workload {
                name: name.to_string(),
                interval: interval,
                bam: bam.to_path_buf(),
            });
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

        let file = File::open(bed.clone()).expect("Failed to open bed file");
        let file = io::BufReader::new(file);
        let last = *extensions.last().unwrap();

        if last == "gz" || last == "gzip" {
            assert_eq!(extensions[extensions.len() - 2], "bed");
            Workload::_from_bed_interval(io::BufReader::new(GzDecoder::new(file)), bam)
        } else {
            assert_eq!(last, "bed");
            Workload::_from_bed_interval(file, bam)
        }
    }

    pub fn len(&self) -> u64 {
        self.interval.range().end - self.interval.range().start
    }
}
