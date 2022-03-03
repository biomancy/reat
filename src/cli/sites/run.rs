use std::collections::HashMap;

use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::cli::shared::args::CoreArgs;
use crate::cli::shared::stranding::Stranding;
use crate::cli::sites::args::SiteArgs;
use crate::core::hooks::engine::REATHooksEngine;
use crate::core::mismatches::site::{REATBatchedSiteMismatches, REATSiteMismatchesBuilder};
use crate::core::rpileup::hts::HTSPileupEngine;
use crate::core::rpileup::ncounters::cnt::{BaseNucCounter, IntervalNucCounter, StrandedNucCounter};
use crate::core::runner::REATRunner;

pub fn run(args: &ArgMatches, mut core: CoreArgs, factory: impl Fn() -> ProgressBar) {
    let args = SiteArgs::new(&mut core, args, &factory);

    // Strander & Hooks don't require any further processing
    let mut strander = args.stranding;
    let hooks: REATHooksEngine<REATBatchedSiteMismatches> = REATHooksEngine::new();

    // Mismatchs builder. Always with prefilter since there are no site-level stats right now
    let builder = REATSiteMismatchesBuilder::new(args.maxwsize, core.refnucpred, args.retain, Some(args.prefilter));

    // Initialize basic counter
    let counter = BaseNucCounter::new(args.maxwsize, core.readfilter, core.trim5, core.trim3);
    let counter = IntervalNucCounter::new(counter);

    let mut saveto = csv::WriterBuilder::new().from_writer(core.saveto);
    match core.stranding {
        Stranding::Unstranded => {
            // Compose strander + pileuper
            let pileuper = HTSPileupEngine::new(core.bamfiles, counter);
            // Launch the processing
            let runner = REATRunner::new(builder, strander, pileuper, hooks);
            shared::run(args.workload, runner, factory(), &mut saveto, HashMap::new()).unwrap();
        }
        Stranding::Stranded(x) => {
            // Remove all stranding algorithm -> they are not required
            strander.clear();
            // Compose strander + pileuper
            let deductor = crate::core::stranding::deduct::DeductStrandByDesign::new(x);
            let pileuper = HTSPileupEngine::new(core.bamfiles, StrandedNucCounter::new(counter, deductor));

            // Launch the processing
            let runner = REATRunner::new(builder, strander, pileuper, hooks);
            shared::run(args.workload, runner, factory(), &mut saveto, HashMap::new()).unwrap();
        }
    };
}

// #[cfg(test)]
// mod test {
//     // use bio_types::genome::Interval;
//     // use mockall::predicate::*;
//     // use mockall::Sequence;
//     //
//     // use crate::core::ncounters::{CountsBufferContent, NucCounterContent, NucCounts};
//     // use crate::core::dna::Nucleotide;
//     // use crate::core::hooks::filters::MockLocusSummaryFilter;
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
