use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::core::filtering::summary::SummaryFilterByMismatches;
use crate::core::workload::ROIWorkload;

use super::args;

pub fn work(pbar: ProgressBar, matches: &ArgMatches) -> (Vec<ROIWorkload>, u32) {
    let roi: &Path = matches.value_of(args::special::ROI).unwrap().as_ref();
    let binsize = matches.value_of(shared::args::core::BINSIZE).unwrap().parse().unwrap();
    pbar.set_message(format!("Parsing BED regions of interest from {}...", roi.display()));
    let workload = ROIWorkload::from_bed(roi, binsize);
    let maxlen = workload.iter().max_by_key(|x| x.len()).map(|x| x.len()).unwrap_or(0);
    pbar.finish_with_message(format!(
        "Will summarize {} ROI editing for regions with max bin size {}",
        workload.len(),
        maxlen
    ));
    (workload, maxlen as u32)
}

pub fn editing_index(pbar: ProgressBar, matches: &ArgMatches) -> Option<BufWriter<File>> {
    matches.value_of(args::stats::EDITING_INDEX).map(|x| {
        pbar.set_message("Parsing EI output path...");
        let ei = PathBuf::from_str(x).unwrap();
        let file = BufWriter::new(File::create(ei.as_path()).unwrap());
        pbar.finish_with_message(format!("Editing indices will be saved to {}", ei.display()));
        file
    })
}
