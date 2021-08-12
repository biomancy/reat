use std::fs::File;

use std::io::BufWriter;
use std::path::Path;
use std::str::FromStr;

use bio_types::genome::Interval;
use clap::{crate_authors, crate_name, crate_version, App, AppSettings, ArgMatches};

use rayon::ThreadPoolBuilder;

use cli::args;

use crate::cli::stranding::Stranding::{Stranded, Unstranded};
use crate::rada::counting::{BaseNucCounter, StrandedCountsBuffer, UnstrandedCountsBuffer};
use crate::rada::filtering::reads::ReadsFilterByQuality;

use crate::rada::filtering::summary::SummaryFilterByMismatches;
use crate::rada::refnuc::RefNucPredByHeurisitc;

use crate::rada::stranding::deduct::DeductStrandByDesign;

use crate::rada::stranding::predict::{NaiveSequentialStrandPredictor, StrandByAtoIEditing, StrandByGenomicFeatures};

use crate::rada::workload::Workload;

mod cli;
mod rada;

pub fn sumfilter(matches: &ArgMatches) -> SummaryFilterByMismatches {
    let (minmismatches, minfreq) = (
        matches.value_of(args::OUT_MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::OUT_MIN_FREQ).unwrap().parse().unwrap(),
    );
    SummaryFilterByMismatches::new(minmismatches, minfreq)
}

pub fn readfilter(matches: &ArgMatches) -> ReadsFilterByQuality {
    let (mapq, phread) = (
        matches.value_of(args::MAPQ).unwrap().parse().unwrap(),
        matches.value_of(args::PHREAD).unwrap().parse().unwrap(),
    );
    ReadsFilterByQuality::new(mapq, phread)
}

pub fn strandpred(matches: &ArgMatches) -> NaiveSequentialStrandPredictor {
    let (minmismatches, minfreq) = (
        matches.value_of(args::STRANDING_MIN_MISMATCHES).unwrap().parse().unwrap(),
        matches.value_of(args::STRANDING_MIN_FREQ).unwrap().parse().unwrap(),
    );
    let strand_by_editing = Some(StrandByAtoIEditing::new(minmismatches, minfreq));
    let strand_by_features = matches
        .value_of(args::ANNOTATION)
        .map(|x| Some(StrandByGenomicFeatures::from_gff3(x.as_ref())))
        .unwrap_or(None);
    NaiveSequentialStrandPredictor::new(strand_by_editing, strand_by_features)
}

pub fn refnucpred(matches: &ArgMatches) -> RefNucPredByHeurisitc {
    let (mincoverage, minfreq) = (
        matches.value_of(args::AUTOREF_MIN_COVERAGE).unwrap().parse().unwrap(),
        matches.value_of(args::AUTOREF_MIN_FREQ).unwrap().parse().unwrap(),
    );
    RefNucPredByHeurisitc::new(mincoverage, minfreq)
}

pub fn workload(bamfiles: &[&Path], matches: &ArgMatches) -> (Vec<Workload>, u32) {
    if let Some(binsize) = matches.value_of(args::BINSIZE) {
        let binsize = binsize.parse().unwrap();
        let workload = Workload::from_binned_hts(bamfiles, binsize);
        (workload, binsize as u32)
    } else {
        let workload = Workload::from_bed_intervals(matches.value_of(args::ROI).unwrap().as_ref());
        let maxlen = workload.iter().max_by_key(|x| x.len()).map(|x| x.len()).unwrap_or(0);
        (workload, maxlen as u32)
    }
}

fn main() {
    let matches = App::new(crate_name!())
        .author(crate_authors!("\n"))
        .version(crate_version!())
        .max_term_width(120)
        .setting(AppSettings::DeriveDisplayOrder)
        .args(cli::args::all())
        .get_matches();

    let threads = matches.value_of(args::THREADS).unwrap().parse().unwrap();
    ThreadPoolBuilder::new().num_threads(threads).build_global().expect("Failed to initialize global thread pool");

    let bamfiles: Vec<&Path> = matches.values_of(args::INPUT).unwrap().map(|x| Path::new(x)).collect();
    let reference = matches.value_of(args::REFERENCE).unwrap().as_ref();

    let (workload, maxsize) = workload(&bamfiles, &matches);

    let readfilter = readfilter(&matches);
    let refnucpred = refnucpred(&matches);
    let strandpred = strandpred(&matches);
    let sumfilter = sumfilter(&matches);

    let dummy = Interval::new("".into(), 1..2);
    let stranding = cli::stranding::Stranding::from_str(matches.value_of(args::STRANDING).unwrap()).unwrap();

    let saveto = matches.value_of(args::SAVETO).unwrap();
    let saveto = BufWriter::new(File::create(saveto).unwrap());

    // TODO: refactor this using a builder-like pattern
    match stranding {
        Unstranded => {
            let counter = BaseNucCounter::new(readfilter, UnstrandedCountsBuffer::new(maxsize), dummy);
            if matches.is_present(args::ROI) {
                cli::resformat::regions(
                    saveto,
                    rada::run::regions(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
                );
            } else {
                cli::resformat::loci(
                    saveto,
                    rada::run::loci(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
                );
            }
        }
        Stranded(design) => {
            let deductor = DeductStrandByDesign::new(design);
            let counter = BaseNucCounter::new(readfilter, StrandedCountsBuffer::new(maxsize, deductor), dummy);
            if matches.is_present(args::ROI) {
                cli::resformat::regions(
                    saveto,
                    rada::run::regions(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
                );
            } else {
                cli::resformat::loci(
                    saveto,
                    rada::run::loci(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
                );
            }
        }
    }
}
