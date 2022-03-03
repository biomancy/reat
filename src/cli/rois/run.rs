use std::collections::HashMap;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::rois::args::ROIArgs;
use crate::cli::shared;
use crate::cli::shared::args::CoreArgs;
use crate::cli::shared::stranding::Stranding;
use crate::core::hooks::engine::REATHooksEngine;
use crate::core::hooks::filters;
use crate::core::hooks::stats::ROIEditingIndex;
use crate::core::mismatches::roi::{REATBatchedROIMismatches, REATROIMismatchesBuilder};
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::cnt::{BaseNucCounter, ROINucCounter, StrandedNucCounter};
use crate::core::runner::REATRunner;

pub fn run(args: &ArgMatches, core: CoreArgs, factory: impl Fn() -> ProgressBar) {
    let args = ROIArgs::new(&core, args, &factory);

    let mut hooks: REATHooksEngine<REATBatchedROIMismatches> = REATHooksEngine::new();
    let builder = match args.ei {
        None => {
            // Always with prefilter since there are no site-level stats right now
            REATROIMismatchesBuilder::new(args.maxwsize, core.refnucpred, args.retain, Some(args.prefilter))
        }
        Some(_) => {
            // Disable prefilter and use a hook instead
            hooks.add_stat(Box::new(ROIEditingIndex::default()));
            let filter: filters::ByMismatches = args.prefilter.into();
            hooks.add_filter(Box::new(filter));
            // Builder without prefiltering
            REATROIMismatchesBuilder::new(args.maxwsize, core.refnucpred, args.retain, None)
        }
    };

    // Initialize basic counter
    let counter = BaseNucCounter::new(args.maxwsize, core.readfilter, core.trim5, core.trim3);
    let counter = ROINucCounter::new(counter);

    let mut saveto = csv::WriterBuilder::new().from_writer(core.saveto);

    let mut strander = args.stranding;
    match core.stranding {
        Stranding::Unstranded => {
            // Compose strander + pileuper
            let pileuper = HTSPileupEngine::new(core.bamfiles, counter);
            // Launch the processing
            let runner = REATRunner::new(builder, strander, pileuper, hooks);
            shared::run(args.workload, runner, factory(), &mut saveto, HashMap::new()).unwrap()
        }
        Stranding::Stranded(x) => {
            // Remove all stranding algorithm -> they are not required
            strander.clear();
            // Compose strander + pileuper
            let deductor = crate::core::stranding::deduct::DeductStrandByDesign::new(x);
            let pileuper = HTSPileupEngine::new(core.bamfiles, StrandedNucCounter::new(counter, deductor));

            // Launch the processing
            let runner = REATRunner::new(builder, strander, pileuper, hooks);
            shared::run(args.workload, runner, factory(), &mut saveto, HashMap::new()).unwrap()
        }
    };

    // debug_assert!(stats.is_empty());

    // resformat::rois(&mut core.saveto, rois);
}

// fn process<
//     Rnr: Runner + Clone + Send,
//     StrandPred: ROIStrandPredictor + Clone + Send,
//     Filter: FilteringEngine,
//     StatBuilder: ROIBasedStat + Clone + Send,
//     OnIteration: Fn(&[ROISummary]) + Sync,
//     OnFinish: FnOnce(&[ROISummary], u32),
// >(
//     workload: Vec<ROIWorkload>,
//     runner: Rnr,
//     stranding: StrandPred,
//     filter: Filter,
//     statbuilder: Option<StatBuilder>,
//     oniter: OnIteration,
//     onfinish: OnFinish,
// ) -> (Vec<ROISummary>, Option<StatBuilder>) {
//     let dostat = statbuilder.is_some();
//     let ctxstore = ThreadCache::new(move || {
//         RefCell::new((runner.clone(), stranding.clone(), filter.clone(), statbuilder.clone()))
//     });
//
//     let results: Vec<ROISummary> = workload
//         .into_par_iter()
//         .map(|w| {
//             let mut ctx = ctxstore.get().borrow_mut();
//             let (rnr, stranding, filter, stat) = ctx.deref_mut();
//             let result = doiter(w, rnr, stranding, filter, stat);
//             oniter(&result);
//             result
//         })
//         .flatten()
//         .collect();
//
//     let ctxstore = ctxstore.dissolve().map(|x| x.into_inner()).collect_vec();
//
//     let mapped = ctxstore.iter().map(|x| x.0.mapped()).sum();
//
//     let stats = if dostat {
//         let buffer = ctxstore.into_iter().map(|x| x.3.unwrap()).collect_vec();
//         Some(StatBuilder::combine(buffer.as_slice()))
//     } else {
//         None
//     };
//
//     onfinish(&results, mapped);
//
//     (results, stats)
// }
//
// fn doiter(
//     w: ROIWorkload,
//     rnr: &mut impl Runner,
//     stranding: &impl ROIStrandPredictor,
//     filter: &impl FilteringEngine,
//     stat: &mut Option<impl ROIBasedStat>,
// ) -> Vec<ROISummary> {
//     // Count nucleotides and infer the reference sequence / strand
//     let (mut unstranded, stranded): (Vec<ROISummary>, Vec<ROISummary>) =
//         rnr.run(w).into_iter().map(|x| x.into()).partition(|x: &ROISummary| x.strand.is_unknown());
//
//     let strands = stranding.batch_predict(&unstranded);
//     for (item, strand) in zip(&mut unstranded, strands) {
//         item.strand = strand;
//     }
//     let mut results = Vec::with_capacity(stranded.len() + unstranded.len());
//     results.extend(unstranded.into_iter());
//     results.extend(stranded.into_iter());
//
//     if let Some(stat) = stat {
//         for x in &results {
//             stat.process(x);
//         }
//     }
//
//     results.retain(|x| filter.is_ok(x));
//     results
// }

// #[cfg(test)]
// mod test {
//     // use bio_types::genome::Interval;
//     // use mockall::predicate::*;
//     // use mockall::Sequence;
//     //
//     // use crate::core::ncounters::{CountsBufferContent, NucCounterContent, NucCounts};
//     // use crate::core::dna::Nucleotide;
//     // use crate::core::hooks::filters::MockIntervalSummaryFilter;
//     // use crate::core::run::runner::MockROIRunCtx;
//     // use crate::core::stats::MockIntervalBasedStat;
//     // use crate::core::stranding::predict::MockIntervalStrandPredictor;
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
//     // fn mock_strandpred(returning: Strand) -> MockIntervalStrandPredictor {
//     //     let mut mock = MockIntervalStrandPredictor::new();
//     //     mock.expect_predict().once().return_const(returning);
//     //     mock
//     // }
//     //
//     // fn mock_filter(returning: bool) -> MockIntervalSummaryFilter {
//     //     let mut mock = MockIntervalSummaryFilter::new();
//     //     mock.expect_is_ok().once().return_const(returning);
//     //     mock
//     // }
//     //
//     // #[test]
//     // fn run_empty() {
//     //     let mut runner = MockROIRunCtx::new();
//     //     let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };
//     //
//     //     // Empty result
//     //     let mut seq = Sequence::new();
//     //     runner.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
//     //     runner.expect_finalize()
//     //         .once()
//     //         .in_sequence(&mut seq)
//     //         .returning(|| Option::<(NucCounterContent, Vec<Nucleotide>)>::None);
//     //     let mut stat = MockIntervalBasedStat::new();
//     //     assert!(super::_run(&w, &mut runner, Some(&mut stat)).is_empty());
//     // }
//     //
//     // #[test]
//     // fn run_stranded() {
//     //     let mut runner = MockROIRunCtx::new();
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
//     //
//     //     runner.expect_finalize()
//     //         .once()
//     //         .in_sequence(&mut seq)
//     //         .return_const(Some((content, vec![Nucleotide::G, Nucleotide::T])));
//     //     let mut stat = MockIntervalBasedStat::new();
//     //     stat.expect_process().times(2).in_sequence(&mut seq).return_const(());
//     //     for filter in [true, false] {
//     //         runner.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
//     //     }
//     //     let result = super::_run(&w, &mut runner, Some(&mut stat));
//     //     assert_eq!(result.len(), 1);
//     //     assert!(!result[0].strand.is_unknown());
//     // }
//     //
//     // #[test]
//     // fn run_unstranded() {
//     //     let mut runner = MockROIRunCtx::new();
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
//     //     runner.expect_strandpred().once().in_sequence(&mut seq).return_const(mock_strandpred(Strand::Forward));
//     //     let mut stat = MockIntervalBasedStat::new();
//     //     stat.expect_process().once().in_sequence(&mut seq).return_const(());
//     //     runner.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(true));
//     //
//     //     let result = super::_run(&w, &mut runner, Some(&mut stat));
//     //     assert_eq!(result.len(), 1);
//     //     assert_eq!(result[0].strand, Strand::Forward);
//     // }
// }
