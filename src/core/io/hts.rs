use std::collections::HashMap;
use std::path::Path;

use bio_types::genome::Interval;
use itertools::{zip, Itertools};
pub use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;

pub fn contigs(hts: &[impl AsRef<Path>]) -> Vec<Interval> {
    let mut contigs = HashMap::new();

    let readers = hts
        .iter()
        .map(|file| {
            let file = file.as_ref();
            IndexedReader::from_path(file).unwrap_or_else(|_| {
                panic!(
                    "Failed to open file {}\n\
                        Possible reasons: BAM file was not indexed (samtools index); you don't have read permissions",
                    file.display()
                )
            })
        })
        .collect_vec();
    let headers = readers.iter().map(|x| x.header()).collect_vec();

    for (file, header) in zip(hts, headers) {
        for tid in 0..header.target_count() {
            let name = String::from_utf8_lossy(header.tid2name(tid));
            let length = header
                .target_len(tid)
                .unwrap_or_else(|| panic!("Failed to parse header for {}", file.as_ref().display()));

            let stored = contigs.entry(name.clone()).or_insert(length);
            assert_eq!(
                *stored, length,
                "BAM headers must contain equivalent contigs, {} two lengths {} != {}",
                name, length, stored
            );
        }
    }

    contigs.into_iter().map(|(name, length)| Interval::new(name.into(), 0..length)).collect()
}
