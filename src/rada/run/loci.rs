use std::path::Path;

use bio_types::genome::{AbstractInterval, Locus};
use bio_types::strand::Strand;
use rust_htslib::bam::record::Record;

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::LocusSummaryFilter;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::stranding::predict::LocusStrandPredictor;
use crate::rada::summary::LocusSummary;
use crate::rada::workload::Workload;

use super::context::ThreadContext;

pub fn run<
    Counter: NucCounter<Record>,
    RefNucPred: RefNucPredictor,
    StrandPred: LocusStrandPredictor,
    Filter: LocusSummaryFilter,
>(
    workload: Vec<Workload>,
    bamfiles: &[&Path],
    reference: &Path,
    counter: Counter,
    refnucpred: RefNucPred,
    strandpred: StrandPred,
    filter: Filter,
) -> Vec<LocusSummary> {
    let mut ctx = ThreadContext::new(bamfiles, reference, counter, refnucpred, strandpred, filter);
    workload
        .into_iter()
        .map(|w| {
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

            // Filter results
            result.retain(|x| ctx.filter.is_ok(x));

            result
        })
        .flatten()
        .collect()
}
