use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use bio_types::genome::Interval;
use clap::ArgMatches;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use itertools::Itertools;
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
    saveto: PathBuf,
}

impl ParsedArgs {
    fn new(args: &ArgMatches, mbar: &MultiProgress, factory: impl Fn() -> ProgressBar + Sync) -> ParsedArgs {
        // TODO: is there any other way to please the borrow checker?
        let mut bamfiles: Option<Vec<PathBuf>> = Default::default();
        let mut reference: Option<PathBuf> = Default::default();
        let mut stranding: Option<Stranding> = Default::default();
        let mut readfilter: Option<ReadsFilterByQuality> = Default::default();
        let mut refnucpred: Option<RefNucPredByHeurisitc> = Default::default();
        let mut strandpred: Option<NaiveSequentialStrandPredictor> = Default::default();
        let mut sumfilter: Option<SummaryFilterByMismatches> = Default::default();
        let mut workload: Option<Vec<Workload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();
        let mut saveto: Option<PathBuf> = Default::default();

        let mut pbars = (0..9).map(|_| factory().with_message("....")).rev().collect_vec();

        // TODO: ProgressBar can be only created inside main thread. Yet, docs claim otherwise... Am I missing something?
        // TODO: Refactor
        let (pbarbams, pbarwork) = (pbars.pop().unwrap(), pbars.pop().unwrap());
        rayon::scope(|s| {
            s.spawn(|_| {
                let bams = parse::bamfiles(pbarbams, args);
                let (work, ms) = parse::workload(pbarwork, &bams, args);
                bamfiles = Some(bams);
                workload = Some(work);
                maxsize = Some(ms);
            });
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| reference = Some(parse::reference(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| stranding = Some(parse::stranding(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| readfilter = Some(parse::readfilter(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| refnucpred = Some(parse::refnucpred(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| strandpred = Some(parse::strandpred(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| sumfilter = Some(parse::sumfilter(pbar, args)));
            let pbar = pbars.pop().unwrap();
            s.spawn(|_| saveto = Some(parse::saveto(pbar, args)));
            mbar.join().expect("Failed to render progress bar");
        });
        ParsedArgs {
            bamfiles: bamfiles.unwrap(),
            reference: reference.unwrap(),
            stranding: stranding.unwrap(),
            readfilter: readfilter.unwrap(),
            refnucpred: refnucpred.unwrap(),
            strandpred: strandpred.unwrap(),
            sumfilter: sumfilter.unwrap(),
            workload: workload.unwrap(),
            maxsize: maxsize.unwrap(),
            saveto: saveto.unwrap(),
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

        let args = ParsedArgs::new(matches, &mbar, pbar);
        App { args, matches, mbar }
    }

    fn _run(self, counter: impl NucCounter<Record> + Send + Clone) {
        let style = ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-")
            .on_finish(ProgressFinish::AndLeave);

        let roimode = self.matches.is_present(args::ROI);
        let worksize = if roimode {
            self.args.workload.len()
        } else {
            self.args.workload.len() * self.matches.value_of(args::BINSIZE).unwrap().parse::<usize>().unwrap()
        };
        let pbar = self.mbar.add(ProgressBar::new(worksize as u64).with_style(style));
        pbar.set_draw_delta(100);

        rayon::scope(|s| {
            let args = self.args;

            let saveto = BufWriter::new(File::create(&args.saveto).unwrap());

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
                            |x| pbar.inc(x.len() as u64),
                            || pbar.finish_with_message("FADSAD"),
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
                            |x| pbar.inc(x.len() as u64),
                            || pbar.finish_with_message("FADSAD"),
                        ),
                    );
                });
            }
            self.mbar.join().expect("Failed to render progress bar");
        })
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
