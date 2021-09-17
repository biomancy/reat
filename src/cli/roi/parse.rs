use std::fs::{File, OpenOptions};
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;

use crate::core::workload::ROIWorkload;

use super::args;
use crate::core::stats::EditingIndex;

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
use crate::cli::roi::resformat;
pub fn editing_index(pbar: ProgressBar, matches: &ArgMatches) -> Option<BufWriter<File>> {
    pbar.set_message("Parsing EI output path...");
    match matches.value_of(args::stats::EDITING_INDEX) {
        None => {
            pbar.finish_with_message("Editing index won't be calculated");
            None
        }
        Some(x) => {
            let ei = PathBuf::from_str(x).unwrap();

            let mut options = OpenOptions::new();
            options.write(true);
            let file = match ei.exists() {
                true => BufWriter::new(options.append(true).open(ei.as_path()).unwrap()),
                false => {
                    let mut file = BufWriter::new(options.create(true).open(ei.as_path()).unwrap());
                    resformat::statheader::<EditingIndex, BufWriter<File>>(&mut file);
                    file
                }
            };

            pbar.finish_with_message(format!("Editing index will be saved to {}", ei.display()));
            Some(file)
        }
    }
}
