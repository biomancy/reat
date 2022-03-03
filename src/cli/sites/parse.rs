use std::convert::TryInto;
use std::path::Path;

use bio_types::genome::AbstractInterval;
use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::cli::sites::args::output_filtering::FORCE_LIST;
use crate::core::io;
use crate::core::io::bed;
use crate::core::io::bed::BedRecord;
use crate::core::mismatches::prefilters::retain::RetainSitesFromIntervals;
use crate::core::workload::SiteWorkload;

pub fn work(
    pbar: ProgressBar,
    bamfiles: &[impl AsRef<Path>],
    exclude: Option<Vec<BedRecord>>,
    matches: &ArgMatches,
) -> (Vec<SiteWorkload>, usize) {
    let binsize = matches.value_of(shared::args::core::BINSIZE).unwrap().parse().unwrap();
    pbar.set_message(format!("Splitting the genome into {}bp bins...", binsize));

    let contigs = io::hts::chromosomes(bamfiles);
    let workload = SiteWorkload::from_intervals(contigs, binsize, exclude);
    debug_assert!(!workload.is_empty());

    let maxsize = workload.iter().map(|x| x.range().end - x.range().start).max().unwrap();
    pbar.finish_with_message(format!(
        "Will summarize editing for {} genome bins with max bin size {}",
        workload.len(),
        binsize
    ));
    (workload, maxsize.try_into().unwrap())
}

pub fn retain(pbar: ProgressBar, matches: &ArgMatches) -> Option<RetainSitesFromIntervals> {
    pbar.set_message("Parsing the \"force\" BED file...");

    let bedrecords = matches.value_of(FORCE_LIST).map(|x| bed::parse(Path::new(x)));

    match bedrecords {
        None => {
            pbar.finish_with_message("Forced output is disabled");
            None
        }
        Some(bed) => {
            let loci = bed.iter().map(|x| x.range().end - x.range().start).sum::<u64>();
            pbar.finish_with_message(format!("Output thresholds are disabled for {} sites(force list)", loci));
            Some(RetainSitesFromIntervals::new(bed))
        }
    }
}
