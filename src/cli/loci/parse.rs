use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::loci::args::output_filtering::FORCE_LIST;
use crate::cli::shared;
use crate::core::workload::ROIWorkload;
use std::collections::{HashMap, HashSet};

pub fn work(pbar: ProgressBar, bamfiles: &[impl AsRef<Path>], matches: &ArgMatches) -> (Vec<ROIWorkload>, u32) {
    let binsize = matches.value_of(shared::args::core::BINSIZE).unwrap().parse().unwrap();
    pbar.set_message(format!("Splitting the genome into {}bp bins...", binsize));
    let workload = ROIWorkload::from_hts(bamfiles, binsize);
    pbar.finish_with_message(format!(
        "Will summarize loci editing for {} genome bins with max bin size {}",
        workload.len(),
        binsize
    ));
    (workload, binsize as u32)
}

pub fn forcein(pbar: ProgressBar, matches: &ArgMatches) -> Option<HashMap<String, HashSet<u64>>> {
    pbar.set_message("Parsing the \"force\" BED file...");
    if !matches.is_present(FORCE_LIST) {
        pbar.finish_with_message("Force output of selected regions is disabled");
        return None;
    }

    let mut loci = 0;
    let mut index = HashMap::new();

    let file: String = matches.value_of(FORCE_LIST).unwrap().parse().unwrap();
    let reader = File::open(file).expect("Failed to open bed file");
    let mut reader = io::BufReader::new(reader);

    let mut buf = String::new();
    while reader.read_line(&mut buf).expect("Failed to read BED file") != 0 {
        let line = buf.trim_end();
        if line.is_empty() {
            buf.clear();
            continue;
        }
        let split: Vec<&str> = line.split('\t').take(4).collect();
        assert!(split.len() >= 3);

        let start = split[1].parse::<u64>().expect("Failed to reads string start");
        let end = split[2].parse::<u64>().expect("Failed to reads string start");
        let chrom = split[0].to_owned();

        assert!(end > start);
        for ind in start..end {
            index.entry(chrom.clone()).or_insert(HashSet::new()).insert(ind);
        }
        loci += end - start;
        buf.clear();
    }

    pbar.finish_with_message(format!("Output thresholds are disabled for {} loci(force list)", loci));
    Some(index)
}
