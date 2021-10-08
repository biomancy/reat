use std::cell::RefCell;
use std::ops::DerefMut;

use bio_types::genome::{AbstractInterval, Locus};
use clap::ArgMatches;
use indicatif::ProgressBar;
use itertools::zip;
use rayon::prelude::*;

use crate::cli::loci::args::LociArgs;
use crate::cli::loci::resformat;
use crate::cli::shared;
use crate::cli::shared::args::CoreArgs;
use crate::cli::shared::stranding::Stranding;
use crate::cli::shared::thread_cache::ThreadCache;
use crate::core::counting::buffers::FlatBuffer;
use crate::core::counting::nuccounter::{BaseNucCounter, StrandedNucCounter};
use crate::core::filtering::summary::LocusSummaryFilter;
use crate::core::runner::BaseRunner;
use crate::core::runner::Runner;
use crate::core::stranding::deduct::DeductStrandByDesign;
use crate::core::stranding::predict::LocusStrandPredictor;
use crate::core::summary::LocusSummary;
use crate::core::workload::ROIWorkload;

pub fn run(args: &ArgMatches, mut core: CoreArgs, factory: impl Fn() -> ProgressBar) {
    let args = LociArgs::new(&core, args, &factory);

    // Callbacks to track progress
    let pbar = factory();
    pbar.set_style(shared::style::run::running());
    pbar.set_draw_delta((core.threads * 10) as u64);
    pbar.set_length(args.workload.len() as u64);

    let oniter = |_: &[LocusSummary]| pbar.inc(1);
    let onfinish = |intervals: &[LocusSummary], reads: u32| {
        pbar.set_style(shared::style::run::finished());
        pbar.finish_with_message(format!("Finished with {} loci, total processed reads: {}", intervals.len(), reads))
    };

    let summary = match core.stranding {
        Stranding::Unstranded => {
            let counter =
                BaseNucCounter::new(core.readfilter, FlatBuffer::new(args.maxwsize as usize), core.trim5, core.trim3);
            let runner = BaseRunner::new(core.bamfiles, core.reference, counter, core.refnucpred);
            process(args.workload, runner, args.strandpred, args.outfilter, oniter, onfinish)
        }
        Stranding::Stranded(design) => {
            let deductor = DeductStrandByDesign::new(design);
            let forward = BaseNucCounter::new(
                core.readfilter.clone(),
                FlatBuffer::new(args.maxwsize as usize),
                core.trim5,
                core.trim3,
            );
            let reverse =
                BaseNucCounter::new(core.readfilter, FlatBuffer::new(args.maxwsize as usize), core.trim5, core.trim3);
            let counter = StrandedNucCounter::new(forward, reverse, deductor);
            let runner = BaseRunner::new(core.bamfiles, core.reference, counter, core.refnucpred);
            process(args.workload, runner, args.strandpred, args.outfilter, oniter, onfinish)
        }
    };

    resformat::loci(&mut core.saveto, summary);
}

pub fn process<
    Rnr: Runner + Clone + Send,
    StrandPred: LocusStrandPredictor + Clone + Send,
    Filter: LocusSummaryFilter + Clone + Send,
    OnIteration: Fn(&[LocusSummary]) + Sync,
    OnFinish: FnOnce(&[LocusSummary], u32),
>(
    workload: Vec<ROIWorkload>,
    runner: Rnr,
    strandpred: StrandPred,
    filter: Filter,
    oniter: OnIteration,
    onfinish: OnFinish,
) -> Vec<LocusSummary> {
    let ctxstore = ThreadCache::new(move || RefCell::new((runner.clone(), strandpred.clone(), filter.clone())));
    let results: Vec<LocusSummary> = workload
        .into_par_iter()
        .map(|w| {
            let mut ctx = ctxstore.get().borrow_mut();
            let (rnr, strandpred, filter) = ctx.deref_mut();
            let result = doiter(w, rnr, strandpred, filter);
            oniter(&result);
            result
        })
        .flatten()
        .collect();

    let mapped = ctxstore.dissolve().map(|x| x.into_inner().0.mapped()).sum();
    onfinish(&results, mapped);

    results
}

fn doiter(
    w: ROIWorkload,
    rnr: &mut impl Runner,
    strandpred: &impl LocusStrandPredictor,
    filter: &impl LocusSummaryFilter,
) -> Vec<LocusSummary> {
    let (mut unstranded, stranded): (Vec<LocusSummary>, Vec<LocusSummary>) = rnr
        .run(w)
        .into_iter()
        .filter(|x| x.cnts.coverage > 0)
        .map(|x| {
            let startpos = x.interval.range().start;
            let contig = x.interval.contig();
            zip(x.reference, x.cnts.nuc).enumerate().map(move |(offset, (nuc, cnt))| {
                LocusSummary::new(Locus::new(contig.into(), startpos + offset as u64), x.strand, *nuc, *cnt)
            })
        })
        .flatten()
        .partition(|x: &LocusSummary| x.strand.is_unknown());

    let strands = strandpred.batch_predict(&unstranded);
    for (item, strand) in zip(&mut unstranded, strands) {
        item.strand = strand;
    }

    unstranded.into_iter().chain(stranded.into_iter()).filter(|x| filter.is_ok(x)).collect()
}

// #[cfg(test)]
// mod test {
//     // use bio_types::genome::Interval;
//     // use mockall::predicate::*;
//     // use mockall::Sequence;
//     //
//     // use crate::core::counting::{CountsBufferContent, NucCounterContent, NucCounts};
//     // use crate::core::dna::Nucleotide;
//     // use crate::core::filtering::summary::MockLocusSummaryFilter;
//     // use crate::core::run::runner::MockLociRunCtx;
//     // use crate::core::stranding::predict::MockLocusStrandPredictor;
//     //
//     // use super::*;
//     //
//     // const FORWARD: &'static [NucCounts] =
//     //     &[NucCounts { A: 1, C: 2, G: 12, T: 0 }, NucCounts { A: 0, C: 10, G: 0, T: 0 }];
//     // const REVERSE: &'static [NucCounts] = &[NucCounts { A: 1, C: 0, G: 0, T: 0 }, NucCounts { A: 0, C: 0, G: 0, T: 1 }];
//     // const UNSTRANDED: &'static [NucCounts] = &[
//     //     NucCounts { A: 1, C: 0, G: 0, T: 0 },
//     //     NucCounts { A: 0, C: 1, G: 0, T: 0 },
//     //     NucCounts { A: 0, C: 0, G: 2, T: 0 },
//     //     NucCounts { A: 0, C: 0, G: 2, T: 3 },
//     // ];
//     //
//     // fn mock_strandpred(returning: Strand) -> MockLocusStrandPredictor {
//     //     let mut mock = MockLocusStrandPredictor::new();
//     //     mock.expect_predict().once().return_const(returning);
//     //     mock
//     // }
//     //
//     // fn mock_filter(returning: bool) -> MockLocusSummaryFilter {
//     //     let mut mock = MockLocusSummaryFilter::new();
//     //     mock.expect_is_ok().once().return_const(returning);
//     //     mock
//     // }
//     //
//     // #[test]
//     // fn run_empty() {
//     //     let mut runner = MockLociRunCtx::new();
//     //     let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };
//     //
//     //     // Empty result
//     //     let mut seq = Sequence::new();
//     //     runner.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
//     //     runner.expect_finalize()
//     //         .once()
//     //         .in_sequence(&mut seq)
//     //         .returning(|| Option::<(NucCounterContent, Vec<Nucleotide>)>::None);
//     //     assert!(super::_run(&w, &mut runner).is_empty());
//     // }
//     //
//     // #[test]
//     // fn run_stranded() {
//     //     let mut runner = MockLociRunCtx::new();
//     //     let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };
//     //
//     //     // Stranded results
//     //     let mut seq = Sequence::new();
//     //     runner.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
//     //
//     //     let content = NucCounterContent {
//     //         interval: w.interval.clone(),
//     //         counts: CountsBufferContent { forward: Some(FORWARD), reverse: Some(REVERSE), unstranded: None },
//     //     };
//     //     runner.expect_finalize()
//     //         .once()
//     //         .in_sequence(&mut seq)
//     //         .return_const(Some((content, vec![Nucleotide::G, Nucleotide::T])));
//     //     for filter in [true, true, true, false] {
//     //         runner.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
//     //     }
//     //     let result = super::_run(&w, &mut runner);
//     //     assert_eq!(result.len(), 3);
//     //     let forward = result.iter().filter(|x| x.strand == Strand::Forward).count();
//     //     let reverse = result.iter().filter(|x| x.strand == Strand::Reverse).count();
//     //     assert!((forward == 1 || forward == 2) && (reverse == 1 || reverse == 2) && forward + reverse == 3);
//     // }
//     //
//     // #[test]
//     // fn run_unstranded() {
//     //     let mut runner = MockLociRunCtx::new();
//     //     let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };
//     //
//     //     // Unstranded results
//     //     let mut seq = Sequence::new();
//     //     runner.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
//     //
//     //     let content = NucCounterContent {
//     //         interval: w.interval.clone(),
//     //         counts: CountsBufferContent { forward: None, reverse: None, unstranded: Some(UNSTRANDED) },
//     //     };
//     //     runner.expect_finalize()
//     //         .once()
//     //         .in_sequence(&mut seq)
//     //         .return_const(Some((content, vec![Nucleotide::G, Nucleotide::T, Nucleotide::A, Nucleotide::C])));
//     //     for (strand, filter) in
//     //         zip([Strand::Forward, Strand::Reverse, Strand::Forward, Strand::Reverse], [false, true, true, true])
//     //     {
//     //         runner.expect_strandpred().once().in_sequence(&mut seq).return_const(mock_strandpred(strand));
//     //         runner.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
//     //     }
//     //     let result = super::_run(&w, &mut runner);
//     //     assert_eq!(result.len(), 3);
//     //     let forward = result.iter().filter(|x| x.strand == Strand::Forward).count();
//     //     let reverse = result.iter().filter(|x| x.strand == Strand::Reverse).count();
//     //     assert!(reverse == 2 && forward == 1);
//     // }
// }
