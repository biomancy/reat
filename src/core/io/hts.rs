use std::path::Path;

use bio_types::genome::Interval;
use itertools::Itertools;
pub use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;

pub fn chromosomes(hts: &[impl AsRef<Path>]) -> Vec<Interval> {
    let hts = hts
        .iter()
        .map(|x| {
            IndexedReader::from_path(x).unwrap_or_else(|_| {
                panic!(
                    "Failed to open file {}\n\
                        Possible reasons: BAM file was not indexed (samtools index); you don't have read permissions",
                    x.as_ref().display()
                )
            })
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
    let mut chrs = Vec::new();
    for tid in 0..header.target_count() {
        let tname = String::from_utf8_lossy(header.tid2name(tid)).to_string();
        let tlen = header.target_len(tid).unwrap();
        chrs.push(Interval::new(tname, 0..tlen));
    }
    chrs
}
