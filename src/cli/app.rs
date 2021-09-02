use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bio_types::genome::Interval;
use clap::ArgMatches;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressFinish, ProgressStyle};
use itertools::Itertools;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::Record;

use crate::cli::stranding::Stranding;
use crate::cli::{args, resformat};
use crate::core;
use crate::core::counting::{BaseNucCounter, NucCounter, StrandedCountsBuffer, UnstrandedCountsBuffer};
use crate::core::filtering::reads::ReadsFilterByQuality;
use crate::core::filtering::summary::SummaryFilterByMismatches;
use crate::core::refnuc::RefNucPredByHeurisitc;
use crate::core::run::workload::Workload;
use crate::core::run::{BaseRunCtx, LociRunCtx, ROIRunCtx};
use crate::core::stats::EditingIndex;
use crate::core::stranding::deduct::DeductStrandByDesign;
use crate::core::stranding::predict::NaiveSequentialStrandPredictor;
use crate::core::summary::{IntervalSummary, LocusSummary};

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
    ei: Option<PathBuf>,
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

        let ei =
            if args.is_present(args::STAT_EDITING_INDEX) { Some(parse::editing_index(factory(), args)) } else { None };

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
            ei,
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

    fn _run_roi(
        workload: Vec<Workload>,
        ei: Option<PathBuf>,
        ctx: impl ROIRunCtx + Send + Clone,
        pbar: ProgressBar,
        saveto: &mut impl Write,
    ) {
        let oniter = |_: &[IntervalSummary]| pbar.inc(1);
        let onfinish = |intervals: &[IntervalSummary], reads: usize| {
            pbar.finish_with_message(format!(
                "Finished with {} regions, total counted reads: {}",
                intervals.len(),
                reads
            ))
        };
        match ei {
            None => {
                resformat::regions(
                    saveto,
                    core::run::regions(workload, ctx, Option::<EditingIndex>::None, oniter, onfinish).0,
                );
            }
            Some(eifile) => {
                let files = ctx.htsfiles().iter().map(|x| x.file_name().unwrap().to_str().unwrap()).join(",");
                let (summary, eivalues) =
                    core::run::regions(workload, ctx, Some(EditingIndex::default()), oniter, onfinish);

                resformat::regions(saveto, summary);

                match eifile.exists() {
                    true => {
                        let eifile = &mut BufWriter::new(OpenOptions::new().append(true).open(&eifile).unwrap());
                        resformat::statistic(&files, eifile, &eivalues.unwrap(), false);
                    }
                    false => {
                        let eifile = &mut BufWriter::new(File::create(&eifile).unwrap());
                        resformat::statistic(&files, eifile, &eivalues.unwrap(), true);
                    }
                }
            }
        }
    }

    fn _run_loci(
        workload: Vec<Workload>,
        ctx: impl LociRunCtx + Send + Clone,
        pbar: ProgressBar,
        saveto: &mut impl Write,
    ) {
        let onfinish = |intervals: &[LocusSummary], reads: usize| {
            pbar.finish_with_message(format!(
                "Finished with {} genomic bins; total counted reads: {}",
                intervals.len(),
                reads
            ))
        };
        resformat::loci(saveto, core::run::loci(workload, ctx, |_| pbar.inc(1), onfinish));
    }

    fn _run(self, counter: impl NucCounter<Record> + Send + Clone) {
        let style = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-")
            .on_finish(ProgressFinish::AndLeave);

        let saveto = &mut BufWriter::new(File::create(&self.args.saveto).unwrap());

        let args = self.args;
        // Borrow checker can't handle a complex args usage. But it works with plain old variables
        let (threads, mbar, workload, ei) = (args.threads, self.mbar, args.workload, args.ei);

        let pbar = mbar.add(ProgressBar::new(workload.len() as u64).with_style(style));
        pbar.set_draw_delta((threads * 10) as u64);

        let ctx =
            BaseRunCtx::new(&args.bamfiles, &args.reference, counter, args.refnucpred, args.strandpred, args.sumfilter);

        let roimode = self.matches.is_present(args::ROI);
        rayon::scope(|s| {
            if roimode {
                s.spawn(|_| {
                    Self::_run_roi(workload, ei, ctx, pbar, saveto);
                });
            } else {
                s.spawn(|_| Self::_run_loci(workload, ctx, pbar, saveto));
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
