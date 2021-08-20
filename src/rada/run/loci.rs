use std::cell::RefCell;
use std::sync::Mutex;

use bio_types::genome::{AbstractInterval, Locus};
use bio_types::strand::Strand;
use rayon::prelude::*;
use thread_local::ThreadLocal;

use crate::rada::filtering::summary::LocusSummaryFilter;
use crate::rada::run::ctx::LociRunCtx;
use crate::rada::run::workload::Workload;
use crate::rada::stranding::predict::LocusStrandPredictor;
use crate::rada::summary::LocusSummary;
use itertools::zip;

pub fn run<
    Context: LociRunCtx + Send + Clone,
    OnIteration: Fn(&[LocusSummary]) + Sync,
    OnFinish: FnOnce(&[LocusSummary]),
>(
    workload: Vec<Workload>,
    ctx: Context,
    oniter: OnIteration,
    onfinish: OnFinish,
) -> Vec<LocusSummary> {
    let ctxstore: ThreadLocal<RefCell<Context>> = ThreadLocal::new();
    let ctxbuilder = Mutex::new(move || RefCell::new(ctx.clone()));

    let result: Vec<LocusSummary> = workload
        .into_par_iter()
        .map(|w| {
            let mut ctx = ctxstore.get_or(|| ctxbuilder.lock().unwrap()()).borrow_mut();

            // Count nucleotides for each locus in the given interval
            ctx.count(&w.interval);
            let result = ctx.finalize();
            if result.is_none() {
                let result = vec![];
                oniter(&result);
                return result;
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
                    let locstrand =
                        if strand.is_unknown() { ctx.strandpred().predict(&locus, refnuc, cnt) } else { strand };
                    let summary = LocusSummary::new(locus, locstrand, *refnuc, *cnt);
                    if ctx.filter().is_ok(&summary) {
                        result.push(summary);
                    }
                }
            }

            oniter(&result);
            result
        })
        .flatten()
        .collect();
    onfinish(&result);
    result
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use crate::rada::counting::MockNucCounter;
//     use crate::rada::filtering::summary::MockLocusSummaryFilter;
//     use crate::rada::read::MockRead;
//     use crate::rada::refnuc::MockRefNucPredictor;
//     use crate::rada::stranding::predict::MockLocusStrandPredictor;
//     use bio_types::genome::Interval;
//     use std::str::FromStr;
//
//     #[test]
//     fn run() {
//         let workload = vec![
//             Workload { interval: Interval::new("1".into(), 1..2), name: "First".into() },
//             Workload { interval: Interval::new("2".into(), 3..4), name: "Second".into() },
//         ];
//         let bamfiles = &vec![PathBuf::from_str("A!").unwrap(), PathBuf::from_str("B!").unwrap()];
//         let reference = Path::new("C!");
//         let mut counter = MockNucCounter::<MockRead>::new();
//         let mut refnucpred = MockRefNucPredictor::new();
//         let mut strandpred = MockLocusStrandPredictor::new();
//         let mut filter = MockLocusSummaryFilter::new();
//
//         super::run(workload, bamfiles, reference, counter, refnucpred, strandpred, filter, |_| {}, |_| {});
//     }
// }
