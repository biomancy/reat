use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use bio_types::genome::{AbstractInterval, Locus};
use bio_types::strand::Strand;
use rayon::prelude::*;
use rust_htslib::bam::record::Record;
use thread_local::ThreadLocal;

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::LocusSummaryFilter;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::run::workload::Workload;
use crate::rada::stranding::predict::LocusStrandPredictor;
use crate::rada::summary::LocusSummary;

use super::context::Context;

#[allow(clippy::too_many_arguments)]
pub fn run<
    Counter: NucCounter<Record> + Send + Clone,
    RefNucPred: RefNucPredictor + Send + Clone,
    StrandPred: LocusStrandPredictor + Send + Clone,
    SumFilter: LocusSummaryFilter + Send + Clone,
    OnIteration: Fn(&[LocusSummary]) + Sync,
    OnFinish: FnOnce(&[LocusSummary]) + Sync,
>(
    workload: Vec<Workload>,
    bamfiles: &[PathBuf],
    reference: &Path,
    counter: Counter,
    refnucpred: RefNucPred,
    strandpred: StrandPred,
    filter: SumFilter,
    hook: OnIteration,
    onfinish: OnFinish,
) -> Vec<LocusSummary> {
    #[allow(clippy::type_complexity)]
    let ctxstore: ThreadLocal<RefCell<Context<Record, Counter, RefNucPred, StrandPred, SumFilter>>> =
        ThreadLocal::new();
    let builder = Mutex::new(move || {
        RefCell::new(Context::new(
            bamfiles,
            reference,
            counter.clone(),
            refnucpred.clone(),
            strandpred.clone(),
            filter.clone(),
        ))
    });
    let result: Vec<LocusSummary> = workload
        .into_par_iter()
        .map(|w| {
            let mut ctx = ctxstore.get_or(|| builder.lock().unwrap()()).borrow_mut();
            // Counting nucleotides occurrence
            ctx.nuccount(&w.interval);
            let content = ctx.nuccontent();
            if content.is_none() {
                let result = vec![];
                hook(&result);
                return result;
            }
            let content = content.unwrap();

            let counts = &content.counts;

            // Fetch reference sequence and predict "real" sequence based on the sequenced nucleotides
            let sequence = ctx.predseq(&ctx.reference(&content.interval), &content.counts);

            // Build summaries
            let total = counts.total_counts();
            let mut result = Vec::with_capacity(total);

            for (strand, counts) in [(Strand::Forward, counts.forward), (Strand::Reverse, counts.reverse)] {
                if let Some(counts) = counts {
                    for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
                        let locus = Locus::new(w.interval.contig().into(), w.interval.range().start + offset as u64);
                        let summary = LocusSummary::new(locus, strand, *nuc, *cnt);
                        if ctx.filter.is_ok(&summary) {
                            result.push(summary);
                        }
                    }
                }
            }

            if let Some(counts) = counts.unstranded {
                for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
                    let locus = Locus::new(w.interval.contig().into(), w.interval.range().start + offset as u64);
                    let locstrand = ctx.strandpred.predict(&locus, nuc, cnt);
                    let summary = LocusSummary::new(locus, locstrand, *nuc, *cnt);
                    if ctx.filter.is_ok(&summary) {
                        result.push(summary);
                    }
                }
            }

            hook(&result);
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
