use std::collections::HashMap;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use bio_types::genome::{AbstractInterval, Locus};
use bio_types::strand::{Same, Strand};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::{IntervalSummaryFilter, LocusSummaryFilter};
use crate::rada::read::AlignedRead;
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor, StrandPredictor};
use crate::rada::summary::{IntervalSummary, LocusSummary};
use crate::rada::workload::Workload;

use super::context::ThreadContext;

pub fn run<
    SummaryFilter: LocusSummaryFilter,
    Counter: NucCounter<Record>,
    RefNucPred: RefNucPredictor,
    StrandPred: LocusStrandPredictor,
>(
    ctx: &mut ThreadContext<Record, SummaryFilter, Counter, RefNucPred, StrandPred>,
    work: Workload,
) -> Vec<LocusSummary> {
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
    let mut result = Vec::with_capacity(counts.total_counts());

    for (strand, counts) in [(Strand::Forward, counts.forward), (Strand::Reverse, counts.reverse)] {
        if let Some(counts) = counts {
            for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
                let locus = Locus::new(work.interval.contig().into(), work.interval.range().start + offset as u64);
                result.push(LocusSummary::new(locus, strand, *nuc, *cnt));
            }
        }
    }

    if let Some(counts) = counts.unstranded {
        for (offset, (cnt, nuc)) in counts.iter().zip(&sequence).enumerate() {
            let locus = Locus::new(work.interval.contig().into(), work.interval.range().start + offset as u64);
            let locstrand = ctx.strandpred.predict(&locus, nuc, cnt);
            result.push(LocusSummary::new(locus, locstrand, *nuc, *cnt));
        }
    }

    // Filter results
    result.retain(|x| ctx.filter.is_ok(x));

    result
}
