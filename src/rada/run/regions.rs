use std::cell::RefCell;

use bio_types::strand::Strand;
use rayon::prelude::*;

use crate::rada::filtering::summary::IntervalSummaryFilter;
use crate::rada::run::workload::Workload;
use crate::rada::stats::IntervalBasedStat;
use crate::rada::stranding::predict::IntervalStrandPredictor;
use crate::rada::summary::IntervalSummary;

use super::ctx::ROIRunCtx;
use super::thread_cache::ThreadCache;

pub fn run<
    Context: ROIRunCtx + Send + Clone,
    StatBuilder: IntervalBasedStat + Send + Clone,
    OnIteration: Fn(&[IntervalSummary]) + Sync,
    OnFinish: FnOnce(&[IntervalSummary]),
>(
    workload: Vec<Workload>,
    context: Context,
    stat: Option<StatBuilder>,
    hook: OnIteration,
    onfinish: OnFinish,
) -> (Vec<IntervalSummary>, Option<StatBuilder>) {
    let ctxstore = ThreadCache::new(move || RefCell::new(context.clone()));

    let dostat = stat.is_some();
    let statstore = ThreadCache::new(move || RefCell::new(stat.clone().unwrap()));

    let results: Vec<IntervalSummary> = workload
        .into_par_iter()
        .map(|w| {
            let mut ctx = ctxstore.get().borrow_mut();
            // Count nucleotides and infer the reference sequence
            ctx.count(&w.interval);
            let result = ctx.finalize();
            if result.is_none() {
                let result = vec![];
                hook(&result);
                return result;
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
                    let mut summary = IntervalSummary::from_counts(
                        content.interval.clone(),
                        w.name.clone(),
                        strand,
                        &sequence,
                        counts,
                    );
                    if summary.mismatches.coverage() == 0 {
                        continue;
                    }
                    if strand.is_unknown() {
                        summary.strand = ctx.strandpred().predict(&summary.interval, &summary.mismatches);
                    }
                    result.push(summary);
                }
            }

            if dostat {
                let mut stat = statstore.get().borrow_mut();
                for x in &result {
                    stat.process(x);
                }
            }
            result.retain(|x| ctx.filter().is_ok(x));
            hook(&result);
            result
        })
        .flatten()
        .collect();
    onfinish(&results);

    let stats: Option<StatBuilder> = if dostat {
        Some(statstore.content().map(|x| x.into_inner()).fold(StatBuilder::default(), |a, b| a + b))
    } else {
        None
    };

    (results, stats)
}
