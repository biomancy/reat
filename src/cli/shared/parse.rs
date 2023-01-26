use std::fs::File;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use bio_types::genome::AbstractInterval;
use clap::ArgMatches;
use indicatif::ProgressBar;
use itertools::Itertools;
use rust_htslib::bam::Record;

use crate::cli::shared::stranding::Stranding;
use crate::core::io::fasta::FastaReader;
use crate::core::io::{bed, vcf};
use crate::core::mismatches::{prefilters, MismatchesVec};
use crate::core::refpred::{AutoRef, RefEngine, VCFCorrectedReference};
use crate::core::rpileup::ncounter::filters;
use crate::core::stranding::deduce::StrandSpecificExperimentDesign;
use crate::core::stranding::predict::algo::{StrandByAtoIEditing, StrandByGenomicAnnotation};
use crate::core::stranding::predict::{REATStrandingEngine, StrandingAlgo};

use super::args;

pub fn readfilter(
    pbar: ProgressBar,
    matches: &ArgMatches,
) -> filters::Sequential<Record, filters::ByQuality, filters::ByFlags> {
    pbar.set_message("Parsing filters filter options...");
    let (mapq, nomapq255, phread) = (
        matches.value_of(args::reads_filtering::MAPQ).unwrap().parse().unwrap(),
        matches.is_present(args::reads_filtering::NO_MAPQ_255),
        matches.value_of(args::reads_filtering::PHREAD).unwrap().parse().unwrap(),
    );
    let byquality = filters::ByQuality::new(mapq, nomapq255, phread);

    let (include, exclude) = (
        matches.value_of(args::reads_filtering::INCLUDE_FLAGS).unwrap().parse().unwrap(),
        matches.value_of(args::reads_filtering::EXCLUDE_FLAGS).unwrap().parse().unwrap(),
    );
    let byflags = filters::ByFlags::new(include, exclude);

    let msg = format!(
        "Reads filter options: require flags {}, disallow flags {}, mapq >= {}, phread >= {}. ",
        byflags.include(),
        byflags.exclude(),
        byquality.mapq(),
        byquality.phread()
    );
    if nomapq255 {
        pbar.finish_with_message(msg + "Mapq = 255 is NOT allowed.");
    } else {
        pbar.finish_with_message(msg + "Mapq = 255 is allowed.");
    }

    filters::Sequential::new(byquality, byflags)
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

pub fn saveto(pbar: ProgressBar, matches: &ArgMatches) -> csv::Writer<File> {
    pbar.set_message("Parsing output path...");
    let result = matches.value_of(args::core::SAVETO).unwrap();
    let file = File::create(result).unwrap();
    let file = csv::WriterBuilder::new().from_writer(file);
    pbar.finish_with_message(format!("Result will be saved to {}", result));
    file
}

pub fn stranding(pbar: ProgressBar, matches: &ArgMatches) -> Stranding {
    pbar.set_message("Parsing stranding parameter...");
    let stranding = Stranding::from_str(matches.value_of(args::core::STRANDING).unwrap()).unwrap();
    let msg = match stranding {
        Stranding::Unstranded => {
            "Unstranded library, transcription strand will be predicted by heuristics"
        }
        Stranding::Stranded(x) => {
            match x {
                StrandSpecificExperimentDesign::Same => { "Single-end stranded library: read strand matches transcription strand" }
                StrandSpecificExperimentDesign::Flip => { "Single-end stranded library: read strand is reverse to the transcription strand" }
                StrandSpecificExperimentDesign::Same1Flip2 => { "Paired-end stranded library: read1 matches transcription strand, read2 is reverse to the transcription strand" }
                StrandSpecificExperimentDesign::Flip1Same2 => { "Paired-end stranded library: read1 is reverse to the transcription strand, read2 matches transcription strand" }
            }
        }
    };
    pbar.finish_with_message(msg);
    stranding
}

pub fn strandpred<T>(pbar: ProgressBar, matches: &ArgMatches) -> REATStrandingEngine<T>
where
    T: MismatchesVec,
    StrandByGenomicAnnotation: StrandingAlgo<T>,
    StrandByAtoIEditing: StrandingAlgo<T>,
{
    pbar.set_draw_delta(10_000);
    pbar.set_message("Parsing strand prediction parameters...");

    let mut engine = REATStrandingEngine::new();

    let stranding = Stranding::from_str(matches.value_of(args::core::STRANDING).unwrap()).unwrap();
    if stranding != Stranding::Unstranded {
        pbar.finish_with_message(format!(
            "Strand prediction is disabled -> working with \"{}\" stranded library",
            stranding
        ));
        return engine;
    }

    // User message
    let mut msg = vec![];
    if let Some(x) = matches.value_of(args::stranding::ANNOTATION) {
        msg.push("by genomic features [exons, genes, extended utrs]".to_owned());
        let extend3utr = matches.value_of(args::stranding::EXTEND_UTR3).unwrap_or("0").parse().unwrap();
        engine.add(Box::new(StrandByGenomicAnnotation::from_gff(x.as_ref(), extend3utr, |_| pbar.inc(1))));
    }

    let (minmismatches, minfreq) = (
        matches.value_of(args::stranding::MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::stranding::MIN_FREQ).unwrap().parse().unwrap(),
    );
    msg.push(format!("by A->I editing[min mismatches={}, min freq={}]", minmismatches, minfreq));
    engine.add(Box::new(StrandByAtoIEditing::new(minmismatches, minfreq)));

    let msg = format!("Strand prediction (by priority): {}", msg.join(", "));
    pbar.finish_with_message(msg);
    engine
}

pub fn refnucpred(pbar: ProgressBar, matches: &ArgMatches, reader: Box<dyn FastaReader>) -> Box<dyn RefEngine> {
    pbar.set_message("Parsing reference prediction parameters...");

    if let Some(file) = matches.value_of(args::autoref::VCF) {
        let file = Path::new(file);
        let snv = vcf::parse(file);

        let heterozygotes: usize = snv.heterozygous.iter().map(|x| x.len()).sum();
        let homozygotes: usize = snv.homozygous.iter().map(|x| x.len()).sum();

        let variants = VCFCorrectedReference::new(snv, reader);
        pbar.finish_with_message(format!(
            "Reference will be adjusted by SNPs(heterozygotes: {heterozygotes}, homozygotes: {homozygotes}) from: {}.",
            file.file_name().unwrap().to_str().unwrap()
        ));
        Box::new(variants)
    } else {
        let (mincoverage, minfreq, hyperedit) = (
            matches.value_of(args::autoref::MIN_COVERAGE).unwrap().parse().unwrap(),
            matches.value_of(args::autoref::MIN_FREQ).unwrap().parse().unwrap(),
            matches.is_present(args::autoref::HYPEREDITING),
        );
        let mut msg = format!(
            "Reference prediction for site with coverage >= {} and most common nucleotide frequency >= {}.",
            mincoverage, minfreq
        );
        if hyperedit {
            msg += " A->G or T->C corrections was disabled (hyper editing mode)."
        }
        let result = AutoRef::new(mincoverage, minfreq, hyperedit, reader);
        pbar.finish_with_message(msg);
        Box::new(result)
    }
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
) -> prefilters::ByMismatches {
    pbar.set_message("Parsing filtering options...");
    let (minmismatches, minfreq, mincov) = (
        matches.value_of(mismatch_key).unwrap().parse().unwrap(),
        matches.value_of(freq_key).unwrap().parse().unwrap(),
        matches.value_of(cov_key).unwrap().parse().unwrap(),
    );
    let result = prefilters::ByMismatches::new(minmismatches, minfreq, mincov);
    pbar.finish_with_message(format!(
        "Filtering options: min coverage >= {}; mismatches min number >= {}, min frequency >= {}",
        result.mincov(),
        result.minmismatches(),
        result.minfreq()
    ));
    result
}

pub fn excluded(pbar: ProgressBar, matches: &ArgMatches) -> Option<Vec<bed::BedRecord>> {
    pbar.set_message("Parsing excluded regions...");

    if let Some(path) = matches.value_of(args::core::EXCLUDE_LIST) {
        let bed = bed::parse(Path::new(path));
        let bases = bed.iter().map(|x| x.interval.range().end - x.interval.range().start).sum::<u64>();
        pbar.finish_with_message(format!("Excluded from the processing: {} regions({} bases)", bed.len(), bases));
        Some(bed)
    } else {
        pbar.finish_with_message("No regions will be excluded from the processing");
        None
    }
}
