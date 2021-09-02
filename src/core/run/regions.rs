use std::cell::RefCell;
use std::ops::DerefMut;

use bio_types::strand::Strand;
use rayon::prelude::*;

use crate::core::filtering::summary::IntervalSummaryFilter;
use crate::core::run::workload::Workload;
use crate::core::stats::IntervalBasedStat;
use crate::core::stranding::predict::IntervalStrandPredictor;
use crate::core::summary::ROISummary;

use super::ctx::ROIRunCtx;
use super::thread_cache::ThreadCache;

pub fn run<
    Context: ROIRunCtx + Send + Clone,
    StatBuilder: IntervalBasedStat + Send + Clone,
    OnIteration: Fn(&[ROISummary]) + Sync,
    OnFinish: FnOnce(&[ROISummary], u32),
>(
    workload: Vec<Workload>,
    context: Context,
    stat: Option<StatBuilder>,
    hook: OnIteration,
    onfinish: OnFinish,
) -> (Vec<ROISummary>, Option<StatBuilder>) {
    let ctxstore = ThreadCache::new(move || RefCell::new(context.clone()));

    let statstore =
        if stat.is_some() { Some(ThreadCache::new(move || RefCell::new(stat.clone().unwrap()))) } else { None };

    let results: Vec<ROISummary> = workload
        .into_par_iter()
        .map(|w| {
            let ctx = ctxstore.get();
            let result = if let Some(store) = &statstore {
                _run(&w, ctx.borrow_mut().deref_mut(), Some(store.get().borrow_mut().deref_mut()))
            } else {
                _run(&w, ctx.borrow_mut().deref_mut(), Option::<&mut StatBuilder>::None)
            };
            hook(&result);
            result
        })
        .flatten()
        .collect();

    let stats: Option<StatBuilder> =
        statstore.map(|store| store.content().map(|x| x.into_inner()).fold(StatBuilder::default(), |a, b| a + b));

    let reads_counted = ctxstore.content().map(|x| x.into_inner().reads_counted()).sum();
    onfinish(&results, reads_counted);

    (results, stats)
}

fn _run(w: &Workload, ctx: &mut impl ROIRunCtx, stat: Option<&mut impl IntervalBasedStat>) -> Vec<ROISummary> {
    // Count nucleotides and infer the reference sequence
    ctx.count(&w.interval);
    let result = ctx.finalize();
    if result.is_none() {
        return vec![];
    }
    let (content, sequence) = result.unwrap();

    // Build summaries
    let mut result = Vec::with_capacity(3);
    for (strand, counts) in [
        (Strand::Forward, content.counts.forward),
        (Strand::Reverse, content.counts.reverse),
        (Strand::Unknown, content.counts.unstranded),
    ] {
        if let Some(counts) = counts {
            let mut summary =
                ROISummary::from_counts(content.interval.clone(), w.name.clone(), strand, &sequence, counts);
            if summary.mismatches.coverage() == 0 {
                continue;
            }
            if strand.is_unknown() {
                summary.strand = ctx.strandpred().predict(&summary.interval, &summary.mismatches);
            }
            result.push(summary);
        }
    }

    if let Some(stat) = stat {
        for x in &result {
            stat.process(x);
        }
    };
    result.retain(|x| ctx.filter().is_ok(x));
    result
}

#[cfg(test)]
mod test {
    use bio_types::genome::Interval;
    use mockall::predicate::*;
    use mockall::Sequence;

    use crate::core::counting::{CountsBufferContent, NucCounterContent, NucCounts};
    use crate::core::dna::Nucleotide;
    use crate::core::filtering::summary::MockIntervalSummaryFilter;
    use crate::core::run::ctx::MockROIRunCtx;
    use crate::core::stats::MockIntervalBasedStat;
    use crate::core::stranding::predict::MockIntervalStrandPredictor;

    use super::*;

    const FORWARD: &'static [NucCounts] =
        &[NucCounts { A: 1, C: 2, G: 12, T: 0 }, NucCounts { A: 0, C: 10, G: 0, T: 0 }];
    const REVERSE: &'static [NucCounts] = &[NucCounts { A: 1, C: 0, G: 0, T: 0 }, NucCounts { A: 0, C: 0, G: 0, T: 1 }];
    const UNSTRANDED: &'static [NucCounts] = &[
        NucCounts { A: 1, C: 0, G: 0, T: 0 },
        NucCounts { A: 0, C: 1, G: 0, T: 0 },
        NucCounts { A: 0, C: 0, G: 2, T: 0 },
        NucCounts { A: 0, C: 0, G: 2, T: 3 },
    ];

    fn mock_strandpred(returning: Strand) -> MockIntervalStrandPredictor {
        let mut mock = MockIntervalStrandPredictor::new();
        mock.expect_predict().once().return_const(returning);
        mock
    }

    fn mock_filter(returning: bool) -> MockIntervalSummaryFilter {
        let mut mock = MockIntervalSummaryFilter::new();
        mock.expect_is_ok().once().return_const(returning);
        mock
    }

    #[test]
    fn run_empty() {
        let mut ctx = MockROIRunCtx::new();
        let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };

        // Empty result
        let mut seq = Sequence::new();
        ctx.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
        ctx.expect_finalize()
            .once()
            .in_sequence(&mut seq)
            .returning(|| Option::<(NucCounterContent, Vec<Nucleotide>)>::None);
        let mut stat = MockIntervalBasedStat::new();
        assert!(super::_run(&w, &mut ctx, Some(&mut stat)).is_empty());
    }

    #[test]
    fn run_stranded() {
        let mut ctx = MockROIRunCtx::new();
        let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };

        // Stranded results
        let mut seq = Sequence::new();
        ctx.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());

        let content = NucCounterContent {
            interval: w.interval.clone(),
            counts: CountsBufferContent { forward: Some(FORWARD), reverse: Some(REVERSE), unstranded: None },
        };

        ctx.expect_finalize()
            .once()
            .in_sequence(&mut seq)
            .return_const(Some((content, vec![Nucleotide::G, Nucleotide::T])));
        let mut stat = MockIntervalBasedStat::new();
        stat.expect_process().times(2).in_sequence(&mut seq).return_const(());
        for filter in [true, false] {
            ctx.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
        }
        let result = super::_run(&w, &mut ctx, Some(&mut stat));
        assert_eq!(result.len(), 1);
        assert!(!result[0].strand.is_unknown());
    }

    #[test]
    fn run_unstranded() {
        let mut ctx = MockROIRunCtx::new();
        let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };

        // Unstranded results
        let mut seq = Sequence::new();
        ctx.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());

        let content = NucCounterContent {
            interval: w.interval.clone(),
            counts: CountsBufferContent { forward: None, reverse: None, unstranded: Some(UNSTRANDED) },
        };
        ctx.expect_finalize()
            .once()
            .in_sequence(&mut seq)
            .return_const(Some((content, vec![Nucleotide::G, Nucleotide::T, Nucleotide::A, Nucleotide::C])));
        ctx.expect_strandpred().once().in_sequence(&mut seq).return_const(mock_strandpred(Strand::Forward));
        let mut stat = MockIntervalBasedStat::new();
        stat.expect_process().once().in_sequence(&mut seq).return_const(());
        ctx.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(true));

        let result = super::_run(&w, &mut ctx, Some(&mut stat));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].strand, Strand::Forward);
    }
}
