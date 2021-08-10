use std::collections::HashMap;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bio_types::genome::AbstractInterval;
use bio_types::strand::{Same, Strand};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::IntervalSummaryFilter;
use crate::rada::read::AlignedRead;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::stranding::predict::{IntervalStrandPredictor, StrandPredictor};
use crate::rada::summary::IntervalSummary;
use crate::rada::workload::Workload;

use super::context::ThreadContext;

pub fn run<
    SummaryFilter: IntervalSummaryFilter,
    Counter: NucCounter<Record>,
    RefNucPred: RefNucPredictor,
    StrandPred: IntervalStrandPredictor,
>(
    ctx: &mut ThreadContext<Record, SummaryFilter, Counter, RefNucPred, StrandPred>,
    work: Workload,
) -> Vec<IntervalSummary> {
    // Counting nucleotides occurrence
    ctx.nuccount(&work.bam, &work.interval);
    let content = ctx.counter.content();
    let counts = &content.counts;

    if counts.forward.is_none() && counts.reverse.is_none() && counts.unstranded.is_none() {
        return vec![];
    }

    // Fetch reference sequence and predict "real" sequence based on the sequenced nucleotides
    let sequence = ctx.predseq(&ctx.reference(&content.interval), &content.counts);

    // Build summaries
    let mut result = Vec::with_capacity(3);
    for (strand, counts) in
        [(Strand::Forward, counts.forward), (Strand::Reverse, counts.reverse), (Strand::Unknown, counts.unstranded)]
    {
        if let Some(counts) = counts {
            let mut summary =
                IntervalSummary::from_counts(content.interval.clone(), work.name.clone(), strand, &sequence, counts);
            if strand.same(&Strand::Unknown) {
                summary.strand = ctx.strandpred.predict(&summary.interval, &summary.mismatches);
            }
            result.push(summary);
        }
    }

    // Filter results
    result.retain(|x| ctx.filter.is_ok(x));

    result
}
