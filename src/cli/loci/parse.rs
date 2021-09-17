use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::core::filtering::summary::SummaryFilterByMismatches;
use crate::core::workload::ROIWorkload;

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
