use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;
use itertools::Itertools;

use crate::cli::stranding::Stranding;
use crate::rada::filtering::reads::ReadsFilterByQuality;
use crate::rada::filtering::summary::SummaryFilterByMismatches;
use crate::rada::refnuc::RefNucPredByHeurisitc;
use crate::rada::run::workload::Workload;
use crate::rada::stranding::deduct::StrandSpecificExperimentDesign;
use crate::rada::stranding::predict::{NaiveSequentialStrandPredictor, StrandByAtoIEditing, StrandByGenomicFeatures};

use super::args;

pub fn sumfilter(pbar: ProgressBar, matches: &ArgMatches) -> SummaryFilterByMismatches {
    pbar.set_message("Parsing output filter options...");
    let (minmismatches, minfreq) = (
        matches.value_of(args::OUT_MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::OUT_MIN_FREQ).unwrap().parse().unwrap(),
    );
    let result = SummaryFilterByMismatches::new(minmismatches, minfreq);
    pbar.finish_with_message(format!(
        "Output filter options: mismatches min number >= {}, min frequency >= {}",
        result.minmismatches(),
        result.minfreq()
    ));
    result
}

pub fn readfilter(pbar: ProgressBar, matches: &ArgMatches) -> ReadsFilterByQuality {
    pbar.set_message("Parsing reads filter options...");
    let (mapq, phread) = (
        matches.value_of(args::MAPQ).unwrap().parse().unwrap(),
        matches.value_of(args::PHREAD).unwrap().parse().unwrap(),
    );
    let result = ReadsFilterByQuality::new(mapq, phread);
    pbar.finish_with_message(format!("Reads filter options: mapq >= {}, phread >= {}", result.mapq(), result.phread()));
    result
}

pub fn saveto(pbar: ProgressBar, matches: &ArgMatches) -> PathBuf {
    pbar.set_message("Parsing output path...");
    let result = matches.value_of(args::SAVETO).unwrap();
    pbar.finish_with_message(format!("Result will be saved to {}", result));
    result.into()
}

pub fn stranding(pbar: ProgressBar, matches: &ArgMatches) -> Stranding {
    pbar.set_message("Parsing stranding parameter...");
    let stranding = Stranding::from_str(matches.value_of(args::STRANDING).unwrap()).unwrap();
    let msg = match stranding {
        Stranding::Unstranded => {
            "Unstranded library, regions/loci strand will predicted by heuristic"
        }
        Stranding::Stranded(x) => {
            match x {
                StrandSpecificExperimentDesign::Same => {"Single-end stranded library: read strand matches transcription strand"}
                StrandSpecificExperimentDesign::Flip => {"Single-end stranded library: read strand is reverse to the transcription strand"}
                StrandSpecificExperimentDesign::Same1Flip2 => {"Paired-end stranded library: read1 matches transcription strand, read2 is reverse to the transcription strand"}
                StrandSpecificExperimentDesign::Flip1Same2 => {"Paired-end stranded library: read1 is reverse to the transcription strand, read2 matches transcription strand"}
            }
        }
    };
    pbar.finish_with_message(msg);
    stranding
}

pub fn strandpred(pbar: ProgressBar, matches: &ArgMatches) -> NaiveSequentialStrandPredictor {
    pbar.set_message("Parsing strand prediction parameters...");

    let stranding = Stranding::from_str(matches.value_of(args::STRANDING).unwrap()).unwrap();
    if stranding != Stranding::Unstranded {
        pbar.finish_with_message(format!(
            "Strand prediction is disabled -> working with \"{}\" stranded library",
            stranding
        ));
        return NaiveSequentialStrandPredictor::new(None, None);
    }

    let (minmismatches, minfreq) = (
        matches.value_of(args::STRANDING_MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::STRANDING_MIN_FREQ).unwrap().parse().unwrap(),
    );
    let strand_by_editing = Some(StrandByAtoIEditing::new(minmismatches, minfreq));

    pbar.set_draw_delta(10_000);
    let strand_by_features = matches
        .value_of(args::ANNOTATION)
        .map(|x| Some(StrandByGenomicFeatures::from_gff3(x.as_ref(), |_| pbar.inc(1))))
        .unwrap_or(None);
    let result = NaiveSequentialStrandPredictor::new(strand_by_editing, strand_by_features);

    let msg = "Strand prediction";
    match (result.by_features(), result.by_editing()) {
        (None, None) => {
            unreachable!()
        }
        (Some(_), None) => pbar.finish_with_message(format!("{}: by genomic features", msg)),
        (None, Some(x)) => pbar.finish_with_message(format!(
            "{}: by A->I editing[min mismatches={}, min freq={}]",
            msg,
            x.min_mismatches(),
            x.min_freq()
        )),
        (Some(_), Some(x)) => pbar.finish_with_message(format!(
            "{}: by genomic features IF FAIL by A->I editing[min mismatches={}, min freq={}]",
            msg,
            x.min_mismatches(),
            x.min_freq()
        )),
    };
    result
}

pub fn refnucpred(pbar: ProgressBar, matches: &ArgMatches) -> RefNucPredByHeurisitc {
    pbar.set_message("Parsing reference prediction parameters...");

    let (mincoverage, minfreq) = (
        matches.value_of(args::AUTOREF_MIN_COVERAGE).unwrap().parse().unwrap(),
        matches.value_of(args::AUTOREF_MIN_FREQ).unwrap().parse().unwrap(),
    );
    let result = RefNucPredByHeurisitc::new(mincoverage, minfreq);
    pbar.finish_with_message(format!(
        "Reference prediction for loci with coverage >= {} and most common nucleotide frequency >= {}",
        result.mincoverage(),
        result.freqthr()
    ));
    result
}

pub fn workload(pbar: ProgressBar, bamfiles: &[impl AsRef<Path>], matches: &ArgMatches) -> (Vec<Workload>, u32) {
    if let Some(binsize) = matches.value_of(args::BINSIZE) {
        let binsize = binsize.parse().unwrap();
        pbar.set_message(format!("Splitting the genome into {}bp bins...", binsize));
        let workload = Workload::from_binned_hts(bamfiles, binsize);
        pbar.finish_with_message(format!(
            "Will summarize loci editing for {} genome bins with max bin size {}",
            workload.len(),
            binsize
        ));
        (workload, binsize as u32)
    } else {
        let roi: &Path = matches.value_of(args::ROI).unwrap().as_ref();
        pbar.set_message(format!("Parsing BED regions of interest from {}...", roi.display()));
        let workload = Workload::from_bed_intervals(roi);
        let maxlen = workload.iter().max_by_key(|x| x.len()).map(|x| x.len()).unwrap_or(0);
        pbar.finish_with_message(format!(
            "Will summarize {} ROI editing for regions with max size {}",
            workload.len(),
            maxlen
        ));
        (workload, maxlen as u32)
    }
}

pub fn bamfiles(pbar: ProgressBar, matches: &ArgMatches) -> Vec<PathBuf> {
    pbar.set_message("Parsing paths to the input files...");
    let result: Vec<PathBuf> = matches.values_of(args::INPUT).unwrap().map(|x| x.into()).collect();
    if result.len() == 1 {
        pbar.finish_with_message(format!("Input file path: {}", result[0].display()))
    } else {
        let paths = result.iter().map(|x| x.display()).join(" ");
        pbar.finish_with_message(format!("Input files that will be pooled: {}", paths));
    }
    result
}

pub fn reference(pbar: ProgressBar, matches: &ArgMatches) -> PathBuf {
    pbar.set_message("Parsing path to the reference assembly...");
    let result: PathBuf = matches.value_of(args::REFERENCE).unwrap().into();
    pbar.finish_with_message(format!("Path to the reference assembly: {}", result.display()));
    result
}

pub fn threads(pbar: ProgressBar, matches: &ArgMatches) -> usize {
    pbar.set_message("Parsing number of threads allowed to launch...");
    let result = matches.value_of(args::THREADS).and_then(|x| x.parse().ok()).unwrap();
    pbar.finish_with_message(format!("Using thread pool with at most {} threads", result));
    result
}
