use std::cell::RefCell;
use std::ops::DerefMut;
use std::sync::Mutex;

use bio_types::genome::{AbstractInterval, Locus};
use bio_types::strand::Strand;
use itertools::zip;
use rayon::prelude::*;
use thread_local::ThreadLocal;

use crate::core::filtering::summary::LocusSummaryFilter;
use crate::core::run::ctx::LociRunCtx;
use crate::core::run::workload::Workload;
use crate::core::stranding::predict::LocusStrandPredictor;
use crate::core::summary::LocusSummary;

pub fn run<
    Context: LociRunCtx + Send + Clone,
    OnIteration: Fn(&[LocusSummary]) + Sync,
    OnFinish: FnOnce(&[LocusSummary], usize),
>(
    workload: Vec<Workload>,
    ctx: Context,
    oniter: OnIteration,
    onfinish: OnFinish,
) -> Vec<LocusSummary> {
    let ctxstore: ThreadLocal<RefCell<Context>> = ThreadLocal::new();
    let ctxbuilder = Mutex::new(move || RefCell::new(ctx.clone()));

    let results: Vec<LocusSummary> = workload
        .into_par_iter()
        .map(|w| {
            let ctx = ctxstore.get_or(|| ctxbuilder.lock().unwrap()());
            let result = _run(&w, ctx.borrow_mut().deref_mut());
            oniter(&result);
            result
        })
        .flatten()
        .collect();

    let reads_counted = ctxstore.into_iter().map(|x| x.into_inner().reads_counted()).sum();
    onfinish(&results, reads_counted);

    results
}

fn _run(w: &Workload, ctx: &mut impl LociRunCtx) -> Vec<LocusSummary> {
    // Count nucleotides for each locus in the given interval
    ctx.count(&w.interval);
    let result = ctx.finalize();
    if result.is_none() {
        return vec![];
    }
    let (content, sequence) = result.unwrap();

    // Build loci summaries
    let total = content.counts.total_counts();
    let mut result = Vec::with_capacity(total);

    let (contig, workstart) = (w.interval.contig(), w.interval.range().start);
    for (strand, counts) in [
        (Strand::Forward, content.counts.forward),
        (Strand::Reverse, content.counts.reverse),
        (Strand::Unknown, content.counts.unstranded),
    ] {
        if counts.is_none() {
            continue;
        }
        let counts = counts.unwrap();

        for (offset, (cnt, refnuc)) in zip(counts, &sequence).enumerate() {
            let locus = Locus::new(contig.into(), workstart + offset as u64);
            let locstrand = if strand.is_unknown() { ctx.strandpred().predict(&locus, refnuc, cnt) } else { strand };
            let summary = LocusSummary::new(locus, locstrand, *refnuc, *cnt);
            if ctx.filter().is_ok(&summary) {
                result.push(summary);
            }
        }
    }
    result
}

#[cfg(test)]
mod test {
    use bio_types::genome::Interval;
    use mockall::predicate::*;
    use mockall::Sequence;

    use crate::core::counting::{CountsBufferContent, NucCounterContent, NucCounts};
    use crate::core::dna::Nucleotide;
    use crate::core::filtering::summary::MockLocusSummaryFilter;
    use crate::core::run::ctx::MockLociRunCtx;
    use crate::core::stranding::predict::MockLocusStrandPredictor;

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

    fn mock_strandpred(returning: Strand) -> MockLocusStrandPredictor {
        let mut mock = MockLocusStrandPredictor::new();
        mock.expect_predict().once().return_const(returning);
        mock
    }

    fn mock_filter(returning: bool) -> MockLocusSummaryFilter {
        let mut mock = MockLocusSummaryFilter::new();
        mock.expect_is_ok().once().return_const(returning);
        mock
    }

    #[test]
    fn run_empty() {
        let mut ctx = MockLociRunCtx::new();
        let w = Workload { interval: Interval::new("chr".into(), 1..2), name: "".into() };

        // Empty result
        let mut seq = Sequence::new();
        ctx.expect_count().once().with(eq(w.interval.clone())).in_sequence(&mut seq).return_const(());
        ctx.expect_finalize()
            .once()
            .in_sequence(&mut seq)
            .returning(|| Option::<(NucCounterContent, Vec<Nucleotide>)>::None);
        assert!(super::_run(&w, &mut ctx).is_empty());
    }

    #[test]
    fn run_stranded() {
        let mut ctx = MockLociRunCtx::new();
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
        for filter in [true, true, true, false] {
            ctx.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
        }
        let result = super::_run(&w, &mut ctx);
        assert_eq!(result.len(), 3);
        let forward = result.iter().filter(|x| x.strand == Strand::Forward).count();
        let reverse = result.iter().filter(|x| x.strand == Strand::Reverse).count();
        assert!((forward == 1 || forward == 2) && (reverse == 1 || reverse == 2) && forward + reverse == 3);
    }

    #[test]
    fn run_unstranded() {
        let mut ctx = MockLociRunCtx::new();
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
        for (strand, filter) in
            zip([Strand::Forward, Strand::Reverse, Strand::Forward, Strand::Reverse], [false, true, true, true])
        {
            ctx.expect_strandpred().once().in_sequence(&mut seq).return_const(mock_strandpred(strand));
            ctx.expect_filter().once().in_sequence(&mut seq).return_const(mock_filter(filter));
        }
        let result = super::_run(&w, &mut ctx);
        assert_eq!(result.len(), 3);
        let forward = result.iter().filter(|x| x.strand == Strand::Forward).count();
        let reverse = result.iter().filter(|x| x.strand == Strand::Reverse).count();
        assert!(reverse == 2 && forward == 1);
    }
}
