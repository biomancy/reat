use std::fs::File;
use std::io::{stderr, stdout, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;
use std::thread;
use std::thread::sleep;
use std::time::Duration;

use bio_types::genome::Interval;
use clap::{crate_authors, crate_name, crate_version, App, AppSettings, ArgMatches};
use fern;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressFinish, ProgressIterator, ProgressStyle};
use log::info;
use rayon::ThreadPoolBuilder;

use cli::args;

use crate::cli::stranding::Stranding::{Stranded, Unstranded};
use crate::rada::counting::{BaseNucCounter, StrandedCountsBuffer, UnstrandedCountsBuffer};
use crate::rada::filtering::reads::ReadsFilterByQuality;
use crate::rada::filtering::summary::SummaryFilterByMismatches;
use crate::rada::refnuc::{RefNucPredByHeurisitc, RefNucPredictor};
use crate::rada::stranding::deduct::DeductStrandByDesign;
use crate::rada::stranding::predict::{NaiveSequentialStrandPredictor, StrandByAtoIEditing, StrandByGenomicFeatures};

use self::rada::run::workload::Workload;

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

pub fn abvgd(
    mbar: &mut MultiProgress,
    style: &ProgressStyle,
    args: &ArgMatches,
) -> (ReadsFilterByQuality, RefNucPredByHeurisitc, NaiveSequentialStrandPredictor, SummaryFilterByMismatches) {
    let (mut readf, mut refnuc, mut stranding, mut sumf): (
        Option<ReadsFilterByQuality>,
        Option<RefNucPredByHeurisitc>,
        Option<NaiveSequentialStrandPredictor>,
        Option<SummaryFilterByMismatches>,
    ) = (None, None, None, None);
    rayon::scope(|s| {
        let pbar = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        s.spawn(move |_| {
            let pbar = &pbar;
            pbar.set_message("Constructing reads filter...");
            Some(readfilter(args));
            pbar.finish_with_message("Reads filter constructed");
        });

        let pbar = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        s.spawn(move |_| {
            let pbar = &pbar;
            pbar.set_message("Constructing reads filter...");
            Some(refnucpred(args));
            pbar.finish_with_message("Reference nucleotide prediction algorithm constructed");
        });

        let pbar = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        s.spawn(move |_| {
            let pbar = &pbar;
            pbar.set_message("Constructing reads filter...");
            Some(strandpred(args));
            pbar.finish_with_message("Strand prediction algorithm constructed (not used for stranded libraries)");
        });

        let pbar = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        s.spawn(move |_| {
            let pbar = &pbar;
            pbar.set_message("Constructing summary filter");
            Some(sumfilter(args));
            pbar.finish_with_message("Summary filter constructed");
        });
        mbar.join();
    });
    (readf.unwrap(), refnuc.unwrap(), stranding.unwrap(), sumf.unwrap())
}

fn main() {
    let matches = App::new(crate_name!())
        .author(crate_authors!("\n"))
        .version(crate_version!())
        .max_term_width(120)
        .setting(AppSettings::DeriveDisplayOrder)
        .args(cli::args::all())
        .get_matches();
    let mut mbar = MultiProgress::new();
    let style = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner} {msg}")
        .on_finish(ProgressFinish::AndLeave);

    let threads = matches.value_of(args::THREADS).and_then(|x| x.parse().ok()).unwrap();
    ThreadPoolBuilder::new().num_threads(threads).build_global().expect("Failed to initialize global thread pool");

    let bamfiles: Vec<&Path> = matches.values_of(args::INPUT).unwrap().map(|x| Path::new(x)).collect();
    let reference = Path::new(matches.value_of(args::REFERENCE).unwrap());
    let (readfilter, refnucpred, strandpred, sumfilter): (
        ReadsFilterByQuality,
        RefNucPredByHeurisitc,
        NaiveSequentialStrandPredictor,
        SummaryFilterByMismatches,
    ) = abvgd(&mut mbar, &style, &matches);

    rayon::scope(|s| {
        let pb = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        let h1 = s.spawn(move |_| {
            for i in 0..128 {
                pb.set_message(format!("item #{}", i + 1));
                thread::sleep(Duration::from_millis(15));
            }
            pb.finish_with_message("done");
        });

        let pb = mbar.add(ProgressBar::new_spinner().with_style(style.clone()));
        let h2 = s.spawn(move |_| {
            for _ in 0..3 {
                pb.set_position(0);
                for i in 0..128 {
                    pb.set_message(format!("item #{}", i + 1));
                    pb.inc(1);
                    thread::sleep(Duration::from_millis(8));
                }
            }
            pb.finish_with_message("done");
        });
        mbar.join();
    });
    // let (workload, maxsize) = workload(&bamfiles, &matches);
    //
    // let dummy = Interval::new("".into(), 1..2);
    // let stranding = cli::stranding::Stranding::from_str(matches.value_of(args::STRANDING).unwrap()).unwrap();
    //
    // let saveto = matches.value_of(args::SAVETO).unwrap();
    // let saveto = BufWriter::new(File::create(saveto).unwrap());

    // // TODO: refactor this using a builder-like pattern
    // match stranding {
    //     Unstranded => {
    //         let counter = BaseNucCounter::new(readfilter, UnstrandedCountsBuffer::new(maxsize), dummy);
    //         if matches.is_present(args::ROI) {
    //             cli::resformat::regions(
    //                 saveto,
    //                 rada::run::regions(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
    //             );
    //         } else {
    //             cli::resformat::loci(
    //                 saveto,
    //                 rada::run::loci(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
    //             );
    //         }
    //     }
    //     Stranded(design) => {
    //         let deductor = DeductStrandByDesign::new(design);
    //         let counter = BaseNucCounter::new(readfilter, StrandedCountsBuffer::new(maxsize, deductor), dummy);
    //         if matches.is_present(args::ROI) {
    //             cli::resformat::regions(
    //                 saveto,
    //                 rada::run::regions(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
    //             );
    //         } else {
    //             cli::resformat::loci(
    //                 saveto,
    //                 rada::run::loci(workload, &bamfiles, reference, counter, refnucpred, strandpred, sumfilter),
    //             );
    //         }
    //     }
    // }
}
