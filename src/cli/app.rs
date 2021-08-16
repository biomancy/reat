use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use bio_types::genome::Interval;
use clap::ArgMatches;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::Record;

use crate::cli::stranding::Stranding;
use crate::cli::{args, resformat};
use crate::rada;
use crate::rada::counting::{BaseNucCounter, NucCounter, StrandedCountsBuffer, UnstrandedCountsBuffer};
use crate::rada::filtering::reads::ReadsFilterByQuality;
use crate::rada::filtering::summary::SummaryFilterByMismatches;
use crate::rada::refnuc::RefNucPredByHeurisitc;
use crate::rada::run::workload::Workload;
use crate::rada::stranding::deduct::DeductStrandByDesign;
use crate::rada::stranding::predict::NaiveSequentialStrandPredictor;

use super::parse;

struct ParsedArgs {
    bamfiles: Vec<PathBuf>,
    reference: PathBuf,
    stranding: Stranding,
    readfilter: ReadsFilterByQuality,
    refnucpred: RefNucPredByHeurisitc,
    strandpred: NaiveSequentialStrandPredictor,
    sumfilter: SummaryFilterByMismatches,
    workload: Vec<Workload>,
    maxsize: u32,
    threads: usize,
    saveto: PathBuf,
}

impl ParsedArgs {
    fn new(
        args: &ArgMatches,
        threads: usize,
        mbar: &MultiProgress,
        factory: impl Fn() -> ProgressBar + Sync,
    ) -> ParsedArgs {
        // Fast parse trivial options in the main thread
        let bamfiles = parse::bamfiles(factory(), args);
        let reference = parse::reference(factory(), args);
        let stranding = parse::stranding(factory(), args);
        let readfilter = parse::readfilter(factory(), args);
        let refnucpred = parse::refnucpred(factory(), args);
        let sumfilter = parse::sumfilter(factory(), args);
        let saveto = parse::saveto(factory(), args);

        // TODO: is there any other way to please the borrow checker?
        let mut strandpred: Option<NaiveSequentialStrandPredictor> = Default::default();
        let mut workload: Option<Vec<Workload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();

        let (pbar_workload, pbar_stranding) = (factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::workload(pbar_workload, &bamfiles, args);
                workload = Some(w);
                maxsize = Some(m)
            });
            s.spawn(|_| strandpred = Some(parse::strandpred(pbar_stranding, args)));

            if threads != 1 {
                mbar.join().expect("Failed to render progress bar");
            }
        });
        if threads == 1 {
            mbar.join().expect("Failed to render progress bar");
        }

        ParsedArgs {
            bamfiles,
            reference,
            stranding,
            readfilter,
            refnucpred,
            strandpred: strandpred.unwrap(),
            sumfilter,
            workload: workload.unwrap(),
            maxsize: maxsize.unwrap(),
            threads,
            saveto,
        }
    }
}

pub struct App<'a> {
    args: ParsedArgs,
    matches: &'a ArgMatches,
    mbar: MultiProgress,
}

impl<'a> App<'a> {
    pub fn new(matches: &'a ArgMatches) -> App {
        let style = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {spinner} {msg}")
            .tick_strings(&["▹▹▹▹▹", "▸▹▹▹▹", "▹▸▹▹▹", "▹▹▸▹▹", "▹▹▹▸▹", "▹▹▹▹▸", "▪▪▪▪▪"])
            .on_finish(ProgressFinish::AndLeave);

        let mbar = MultiProgress::new();
        let pbar = || mbar.add(ProgressBar::new_spinner().with_style(style.clone()));

        let threads = parse::threads(pbar(), matches);
        ThreadPoolBuilder::new().num_threads(threads).build_global().expect("Failed to initialize thread pool");

        let args = ParsedArgs::new(matches, threads, &mbar, pbar);
        App { args, matches, mbar }
    }

    fn _run(self, counter: impl NucCounter<Record> + Send + Clone) {
        let style = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-")
            .on_finish(ProgressFinish::AndLeave);

        let saveto = &mut BufWriter::new(File::create(&self.args.saveto).unwrap());

        let roimode = self.matches.is_present(args::ROI);
        let worksize = self.args.workload.len();
        let pbar = self.mbar.add(ProgressBar::new(worksize as u64).with_style(style));
        pbar.set_draw_delta((self.args.threads * 10) as u64);

        let threads = self.args.threads;
        let mbar = self.mbar;
        let args = self.args;
        rayon::scope(|s| {
            if roimode {
                s.spawn(|_| {
                    resformat::regions(
                        saveto,
                        rada::run::regions(
                            args.workload,
                            &args.bamfiles,
                            &args.reference,
                            counter,
                            args.refnucpred,
                            args.strandpred,
                            args.sumfilter,
                            |_| pbar.inc(1),
                            |_| pbar.finish_with_message("Finished"),
                        ),
                    );
                });
            } else {
                s.spawn(|_| {
                    resformat::loci(
                        saveto,
                        rada::run::loci(
                            args.workload,
                            &args.bamfiles,
                            &args.reference,
                            counter,
                            args.refnucpred,
                            args.strandpred,
                            args.sumfilter,
                            |_| pbar.inc(1),
                            |_| pbar.finish_with_message("Finished"),
                        ),
                    );
                });
            }
            if threads > 1 {
                mbar.join().expect("Failed to render progress bar");
            }
        });
        if threads == 1 {
            mbar.join().expect("Failed to render progress bar");
        }
    }

    pub fn run(self) {
        let dummy = Interval::new("".into(), 1..2);
        match self.args.stranding {
            Stranding::Unstranded => {
                let counter =
                    BaseNucCounter::new(self.args.readfilter, UnstrandedCountsBuffer::new(self.args.maxsize), dummy);
                self._run(counter);
            }
            Stranding::Stranded(design) => {
                let deductor = DeductStrandByDesign::new(design);
                let counter = BaseNucCounter::new(
                    self.args.readfilter,
                    StrandedCountsBuffer::new(self.args.maxsize, deductor),
                    dummy,
                );
                self._run(counter);
            }
        }
    }
}
