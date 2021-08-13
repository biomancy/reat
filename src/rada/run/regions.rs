use std::cell::RefCell;
use std::path::Path;
use std::sync::Mutex;

use bio_types::strand::{Same, Strand};
use rayon::prelude::*;
use rust_htslib::bam::record::Record;
use thread_local::ThreadLocal;

use indicatif::{MultiProgress, ProgressBar};

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::IntervalSummaryFilter;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::run::workload::Workload;
use crate::rada::stranding::predict::IntervalStrandPredictor;
use crate::rada::summary::IntervalSummary;

use super::context::ThreadContext;

pub fn run<
    Counter: NucCounter<Record> + Send + Clone,
    RefNucPred: RefNucPredictor + Send + Clone,
    StrandPred: IntervalStrandPredictor + Send + Clone,
    Filter: IntervalSummaryFilter + Send + Clone,
>(
    workload: Vec<Workload>,
    bamfiles: &[&Path],
    reference: &Path,
    counter: Counter,
    refnucpred: RefNucPred,
    strandpred: StrandPred,
    filter: Filter,
) -> Vec<IntervalSummary> {
    let ctxstore: ThreadLocal<RefCell<ThreadContext<Record, Counter, RefNucPred, StrandPred, Filter>>> =
        ThreadLocal::new();
    let builder = Mutex::new(move || {
        RefCell::new(ThreadContext::new(
            bamfiles,
            reference,
            counter.clone(),
            refnucpred.clone(),
            strandpred.clone(),
            filter.clone(),
        ))
    });
    let masterbar = MultiProgress::new();
    let mut barstore: ThreadLocal<ProgressBar> = ThreadLocal::new();
    let results = workload
        .into_par_iter()
        .map(|w| {
            let mut ctx = ctxstore.get_or(|| builder.lock().unwrap()()).borrow_mut();
            // Counting nucleotides occurrence
            ctx.nuccount(&w.interval);
            let content = ctx.counter.content();
            let counts = &content.counts;

            if counts.forward.is_none() && counts.reverse.is_none() && counts.unstranded.is_none() {
                return vec![];
            }

            // Fetch reference sequence and predict "real" sequence based on the sequenced nucleotides
            let sequence = ctx.predseq(&ctx.reference(&content.interval), &content.counts);

            // Build summaries
            let mut result = Vec::with_capacity(3);
            for (strand, counts) in [
                (Strand::Forward, counts.forward),
                (Strand::Reverse, counts.reverse),
                (Strand::Unknown, counts.unstranded),
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
                    if strand.same(&Strand::Unknown) {
                        summary.strand = ctx.strandpred.predict(&summary.interval, &summary.mismatches);
                    }
                    result.push(summary);
                }
            }

            // Filter results
            result.retain(|x| ctx.filter.is_ok(x));
            if result.len() != 0 {
                let pbar = barstore.get_or(|| {
                    let pbar = masterbar.add(ProgressBar::new_spinner());
                    // pbar.set_draw_rate(1);
                    pbar
                });
                pbar.inc(result.len() as u64);
            }

            result
        })
        .flatten()
        .collect();

    masterbar.join();
    results
}
