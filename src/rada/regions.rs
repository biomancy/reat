use std::collections::HashMap;
use std::path::{Path, PathBuf};

use bio_types::genome::AbstractInterval;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};

use crate::rada::modules::counting::{CountsBuffer, NucCounter};
use crate::rada::modules::filtering::summary::IntervalSummaryFilter;
use crate::rada::modules::summarization::{IntervalSummarizer, IntervalSummary};
use crate::rada::modules::workload::Workload;

pub struct ThreadContext<
    SummaryFilter: IntervalSummaryFilter,
    Counter: NucCounter,
    Summarizer: IntervalSummarizer,
> {
    htsreaders: HashMap<PathBuf, bam::IndexedReader>,
    filter: SummaryFilter,
    counter: Counter,
    summarizer: Summarizer,
}

impl<SummaryFilter: IntervalSummaryFilter, Counter: NucCounter, Summarizer: IntervalSummarizer>
    ThreadContext<SummaryFilter, Counter, Summarizer>
{
    pub fn new(
        htsfiles: &[&Path],
        filter: SummaryFilter,
        counter: Counter,
        summarizer: Summarizer,
    ) -> Self {
        let mut htsreaders: HashMap<PathBuf, bam::IndexedReader> = HashMap::new();
        for &hts in htsfiles {
            let reader = bam::IndexedReader::from_path(&hts)
                .expect(&format!("Failed to open file {}", hts.display()));
            htsreaders.insert(PathBuf::from(hts), reader);
        }
        ThreadContext {
            htsreaders,
            filter,
            counter,
            summarizer,
        }
    }
}

pub fn run<
    SummaryFilter: IntervalSummaryFilter,
    Counter: NucCounter,
    Summarizer: IntervalSummarizer,
>(
    ctx: &mut ThreadContext<SummaryFilter, Counter, Summarizer>,
    work: Workload,
) -> Vec<IntervalSummary> {
    // 1. reset the counting
    ctx.counter.reset(&work.interval);

    // 2. do counting
    let mut htsreader = ctx
        .htsreaders
        .get_mut(&work.bam)
        .expect("Bug: thread failed to initialize all BAM readers");

    let err = htsreader.fetch((
        work.interval.contig(),
        work.interval.range().start,
        work.interval.range().end,
    ));
    err.expect(&format!(
        "Failed to fetch records for {}:{}-{} (HTS file corrupted?)",
        work.interval.contig(),
        work.interval.range().start,
        work.interval.range().end
    ));

    let mut record = Record::new();
    while let Some(r) = htsreader.read(&mut record) {
        r.expect("Failed to parse record information (HTS file corrupted?)");
        ctx.counter.process(&mut record);
    }

    // 3. summarize stuff
    let mut summaries = ctx.summarizer.summarize(&work.name, ctx.counter.counts());

    // 4. filter it
    summaries.retain(|x| ctx.filter.is_ok(x));

    summaries
}
