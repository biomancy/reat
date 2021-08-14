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

use super::context::ThreadContext;

#[allow(clippy::too_many_arguments)]
pub fn run<
    Counter: NucCounter<Record> + Send + Clone,
    RefNucPred: RefNucPredictor + Send + Clone,
    StrandPred: LocusStrandPredictor + Send + Clone,
    SumFilter: LocusSummaryFilter + Send + Clone,
    Hook: Fn(&[LocusSummary]) + Sync,
    OnFinish: FnOnce(),
>(
    workload: Vec<Workload>,
    bamfiles: &[PathBuf],
    reference: &Path,
    counter: Counter,
    refnucpred: RefNucPred,
    strandpred: StrandPred,
    filter: SumFilter,
    hook: Hook,
    onfinish: OnFinish,
) -> Vec<LocusSummary> {
    #[allow(clippy::type_complexity)]
    let ctxstore: ThreadLocal<RefCell<ThreadContext<Record, Counter, RefNucPred, StrandPred, SumFilter>>> =
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
    let result = workload
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
            let mut result = Vec::with_capacity(counts.total_counts());

            for (strand, counts) in [(Strand::Forward, counts.forward), (Strand::Reverse, counts.reverse)] {
                if let Some(counts) = counts {
                    for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
                        let locus = Locus::new(w.interval.contig().into(), w.interval.range().start + offset as u64);
                        result.push(LocusSummary::new(locus, strand, *nuc, *cnt));
                    }
                }
            }

            if let Some(counts) = counts.unstranded {
                for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
                    let locus = Locus::new(w.interval.contig().into(), w.interval.range().start + offset as u64);
                    let locstrand = ctx.strandpred.predict(&locus, nuc, cnt);
                    result.push(LocusSummary::new(locus, locstrand, *nuc, *cnt));
                }
            }

            hook(&result);

            // Filter results
            result.retain(|x| ctx.filter.is_ok(x));
            result
        })
        .flatten()
        .collect();
    onfinish();
    result
}
