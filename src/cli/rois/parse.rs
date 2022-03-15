use std::fs;
use std::fs::{File};

use std::path::{Path, PathBuf};

use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::rois::args::output_filtering::FORCE_LIST;
use crate::cli::shared;
use crate::core::io;
use crate::core::io::bed;
use crate::core::io::bed::BedRecord;
use crate::core::mismatches::prefilters::retain::RetainROIFromList;
use crate::core::workload::ROIWorkload;

use super::args;

pub fn work(pbar: ProgressBar, matches: &ArgMatches, exclude: Option<Vec<BedRecord>>) -> (Vec<ROIWorkload>, usize) {
    let roi: &Path = matches.value_of(args::special::ROI).unwrap().as_ref();
    let binsize = matches.value_of(shared::args::core::BINSIZE).unwrap().parse().unwrap();
    pbar.set_message(format!("Parsing BED regions of interest from {}...", roi.display()));

    let roi = io::bed::parse(roi);
    let workload = ROIWorkload::from_bed(roi, binsize, exclude);
    let maxlen = workload.iter().max_by_key(|x| x.len()).map(|x| x.len()).unwrap_or(0);
    pbar.finish_with_message(format!(
        "Will summarize {} ROI editing for regions with max bin size {}",
        workload.len(),
        maxlen
    ));
    (workload, maxlen)
}

pub fn editing_index(pbar: ProgressBar, matches: &ArgMatches) -> Option<csv::Writer<File>> {
    pbar.set_message("Parsing EI output path...");
    match matches.value_of(args::stats::EDITING_INDEX) {
        None => {
            pbar.finish_with_message("Editing index won't be calculated");
            None
        }
        Some(x) => {
            let ei = PathBuf::from_str(x).unwrap();

            let mut builder = csv::WriterBuilder::new();
            let mut options = fs::OpenOptions::new();

            if ei.exists() {
                builder.has_headers(false);
                options.append(true);
            } else {
                options.create(true);
            };

            let stream = options.write(true).open(&ei).expect("Failed to open EI index file for writing.");

            let writer = builder.from_writer(stream);

            pbar.finish_with_message(format!("Editing index will be saved to {}", ei.display()));
            Some(writer)
        }
    }
}

pub fn retain(pbar: ProgressBar, matches: &ArgMatches) -> Option<RetainROIFromList> {
    pbar.set_message("Parsing the \"force\" BED file...");

    let bedrecords = matches.value_of(FORCE_LIST).map(|x| bed::parse(Path::new(x)));

    match bedrecords {
        None => {
            pbar.finish_with_message("Forced output is disabled");
            None
        }
        Some(bed) => {
            pbar.finish_with_message(format!("Output thresholds are disabled for {} ROIs(force list)", bed.len()));
            Some(RetainROIFromList::new(bed))
        }
    }
}
