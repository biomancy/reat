use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::str::FromStr;

use clap::ArgMatches;
use indicatif::ProgressBar;
use itertools::Itertools;
use rust_htslib::bam::Record;

use crate::cli::shared::stranding::Stranding;
use crate::core::filtering::reads::{ReadsFilterByFlags, ReadsFilterByQuality, SequentialReadsFilter};
use crate::core::filtering::summary::SummaryFilterByMismatches;
use crate::core::refnuc::RefNucPredByHeurisitc;
use crate::core::stranding::deduct::StrandSpecificExperimentDesign;
use crate::core::stranding::predict::{SequentialStrandPredictor, StrandByAtoIEditing, StrandByGenomicFeatures};

use super::args;

pub fn readfilter(
    pbar: ProgressBar,
    matches: &ArgMatches,
) -> SequentialReadsFilter<Record, ReadsFilterByQuality, ReadsFilterByFlags> {
    pbar.set_message("Parsing reads filter options...");
    let (mapq, allow_mapq_255, phread) = (
        matches.value_of(args::reads_filtering::MAPQ).unwrap().parse().unwrap(),
        matches.is_present(args::reads_filtering::ALLOW_MAPQ_255),
        matches.value_of(args::reads_filtering::PHREAD).unwrap().parse().unwrap(),
    );
    let byquality = ReadsFilterByQuality::new(mapq, allow_mapq_255, phread);

    let (include, exclude) = (
        matches.value_of(args::reads_filtering::INCLUDE_FLAGS).unwrap().parse().unwrap(),
        matches.value_of(args::reads_filtering::EXCLUDE_FLAGS).unwrap().parse().unwrap(),
    );
    let byflags = ReadsFilterByFlags::new(include, exclude);

    let msg = format!(
        "Reads filter options: require flags {}, disallow flags {}, mapq >= {}, phread >= {}. ",
        byflags.include(),
        byflags.exclude(),
        byquality.mapq(),
        byquality.phread()
    );
    if allow_mapq_255 {
        pbar.finish_with_message(msg + "Mapq = 255 is allowed.");
    } else {
        pbar.finish_with_message(msg + "Mapq = 255 is NOT allowed.");
    }

    SequentialReadsFilter::new(byquality, byflags)
}

pub fn trimming(pbar: ProgressBar, matches: &ArgMatches) -> (u16, u16) {
    pbar.set_message("Parsing trimming options...");
    let (trim5, trim3) = (
        matches.value_of(args::reads_filtering::TRIM5).unwrap().parse().unwrap(),
        matches.value_of(args::reads_filtering::TRIM3).unwrap().parse().unwrap(),
    );

    if trim5 != 0 || trim3 != 0 {
        pbar.finish_with_message(format!("Reads trimming: 5`: {}bp; 3`: {}bp.", trim5, trim3));
    } else {
        pbar.finish_with_message("Reads trimming disabled.");
    }
    (trim5, trim3)
}

pub fn saveto(pbar: ProgressBar, matches: &ArgMatches) -> BufWriter<File> {
    pbar.set_message("Parsing output path...");
    let result = matches.value_of(args::core::SAVETO).unwrap();
    let file = BufWriter::new(File::create(result).unwrap());
    pbar.finish_with_message(format!("Result will be saved to {}", result));
    file
}

pub fn stranding(pbar: ProgressBar, matches: &ArgMatches) -> Stranding {
    pbar.set_message("Parsing stranding parameter...");
    let stranding = Stranding::from_str(matches.value_of(args::core::STRANDING).unwrap()).unwrap();
    let msg = match stranding {
        Stranding::Unstranded => {
            "Unstranded library, regions/loci strand will be predicted by heuristic"
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

pub fn strandpred(pbar: ProgressBar, matches: &ArgMatches) -> SequentialStrandPredictor {
    pbar.set_draw_delta(10_000);
    pbar.set_message("Parsing strand prediction parameters...");

    let stranding = Stranding::from_str(matches.value_of(args::core::STRANDING).unwrap()).unwrap();
    if stranding != Stranding::Unstranded {
        pbar.finish_with_message(format!(
            "Strand prediction is disabled -> working with \"{}\" stranded library",
            stranding
        ));
        return SequentialStrandPredictor::new(None, None);
    }

    let (minmismatches, minfreq) = (
        matches.value_of(args::stranding::MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::stranding::MIN_FREQ).unwrap().parse().unwrap(),
    );
    let strand_by_editing = Some(StrandByAtoIEditing::new(minmismatches, minfreq));

    let strand_by_features = matches
        .value_of(args::stranding::ANNOTATION)
        .map(|x| Some(StrandByGenomicFeatures::from_gff3(x.as_ref(), |_| pbar.inc(1))))
        .unwrap_or(None);
    let result = SequentialStrandPredictor::new(strand_by_editing, strand_by_features);

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

    let (mincoverage, minfreq, hyperedit) = (
        matches.value_of(args::autoref::MIN_COVERAGE).unwrap().parse().unwrap(),
        matches.value_of(args::autoref::MIN_FREQ).unwrap().parse().unwrap(),
        matches.is_present(args::autoref::HYPEREDITING),
    );
    let result = RefNucPredByHeurisitc::new(mincoverage, minfreq, hyperedit);
    let mut msg = format!(
        "Reference prediction for loci with coverage >= {} and most common nucleotide frequency >= {}.",
        result.mincoverage(),
        result.freqthr()
    );
    if hyperedit {
        msg += " A->G or T->C corrections disabled (hyper editing mode)."
    }
    pbar.finish_with_message(msg);
    result
}

pub fn bamfiles(pbar: ProgressBar, matches: &ArgMatches) -> Vec<PathBuf> {
    pbar.set_message("Parsing paths to the input files...");
    let result: Vec<PathBuf> = matches.values_of(args::core::INPUT).unwrap().map(|x| x.into()).collect();
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
    let result: PathBuf = matches.value_of(args::core::REFERENCE).unwrap().into();
    pbar.finish_with_message(format!("Path to the reference assembly: {}", result.display()));
    result
}

pub fn threads(pbar: ProgressBar, matches: &ArgMatches) -> usize {
    pbar.set_message("Parsing number of threads allowed to launch...");
    let result = matches.value_of(args::core::THREADS).and_then(|x| x.parse().ok()).unwrap();
    pbar.finish_with_message(format!(
        "Using thread pool with at most {} threads(+ 1 thread to render progress bar)",
        result
    ));
    result
}

pub fn name(pbar: ProgressBar, matches: &ArgMatches) -> String {
    pbar.set_message("Parsing the run title...");
    let result = matches.value_of(args::core::NAME).and_then(|x| x.parse().ok()).unwrap();
    pbar.finish_with_message(format!("Run title: {}", result));
    result
}

pub fn outfilter(
    pbar: ProgressBar,
    mismatch_key: &str,
    freq_key: &str,
    cov_key: &str,
    matches: &ArgMatches,
) -> SummaryFilterByMismatches {
    pbar.set_message("Parsing filtering options...");
    let (minmismatches, minfreq, mincov) = (
        matches.value_of(mismatch_key).unwrap().parse().unwrap(),
        matches.value_of(freq_key).unwrap().parse().unwrap(),
        matches.value_of(cov_key).unwrap().parse().unwrap(),
    );
    let result = SummaryFilterByMismatches::new(minmismatches, minfreq, mincov);
    pbar.finish_with_message(format!(
        "Filtering options: min coverage >= {}; mismatches min number >= {}, min frequency >= {}",
        result.mincov(),
        result.minmismatches(),
        result.minfreq()
    ));
    result
}
