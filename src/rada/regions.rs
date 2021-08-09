use std::collections::HashMap;
use std::path::{Path, PathBuf};

use bio_types::genome::AbstractInterval;
use rust_htslib::bam;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read;

use crate::rada::counting::NucCounter;
use crate::rada::filtering::summary::IntervalSummaryFilter;
use crate::rada::read::AlignedRead;
use crate::rada::summarization::{IntervalSummarizer, IntervalSummary};
use crate::rada::workload::Workload;
use std::marker::PhantomData;

pub struct ThreadContext<
    R: AlignedRead,
    SummaryFilter: IntervalSummaryFilter,
    Counter: NucCounter<R>,
    Summarizer: IntervalSummarizer,
> {
    htsreaders: HashMap<PathBuf, bam::IndexedReader>,
    filter: SummaryFilter,
    counter: Counter,
    summarizer: Summarizer,
    phantom: PhantomData<R>,
}

impl<R: AlignedRead, SummaryFilter: IntervalSummaryFilter, Counter: NucCounter<R>, Summarizer: IntervalSummarizer>
    ThreadContext<R, SummaryFilter, Counter, Summarizer>
{
    pub fn new(htsfiles: &[&Path], filter: SummaryFilter, counter: Counter, summarizer: Summarizer) -> Self {
        let mut htsreaders: HashMap<PathBuf, bam::IndexedReader> = HashMap::new();
        for &hts in htsfiles {
            let reader = bam::IndexedReader::from_path(&hts).expect(&format!("Failed to open file {}", hts.display()));
            htsreaders.insert(PathBuf::from(hts), reader);
        }
        ThreadContext { htsreaders, filter, counter, summarizer, phantom: Default::default() }
    }
}

pub fn run<SummaryFilter: IntervalSummaryFilter, Counter: NucCounter<Record>, Summarizer: IntervalSummarizer>(
    ctx: &mut ThreadContext<Record, SummaryFilter, Counter, Summarizer>,
    work: Workload,
) -> Vec<IntervalSummary> {
    // 1. reset the counting
    ctx.counter.reset(work.interval.clone());

    // 2. do counting
    let htsreader = ctx.htsreaders.get_mut(&work.bam).expect("Bug: thread failed to initialize all BAM readers");

    let err = htsreader.fetch((work.interval.contig(), work.interval.range().start, work.interval.range().end));
    err.expect(&format!(
        "Failed to fetch reads for {}:{}-{} (HTS file corrupted?)",
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
    let mut summaries = ctx.summarizer.summarize(&work.name, ctx.counter.content());

    // 4. filter it
    summaries.retain(|x| ctx.filter.is_ok(x));

    summaries
}
