use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use bio_types::genome::Interval;
use clap::{App, Arg, ArgSettings};
use rayon;
use rayon::prelude::*;
use rust_htslib::faidx;

use crate::rada::modules::counting::{BaseNucCounter, NucCounter, UnstrandedCountsBuffer};
use crate::rada::modules::filtering::records::RecordsFilterByQuality;
use crate::rada::modules::filtering::summary::{IntervalSummaryFilter, SummaryFilterByEditing};
use crate::rada::modules::refnuc::RefNucPredByHeurisitc;
use crate::rada::modules::stranding::predict::{
    DynamicSequentialStrandPredictor, IntervalStrandPredictor, StrandByAtoIEditing,
    StrandByGenomicFeatures,
};
use crate::rada::modules::summarization::{BaseSummarizer, IntervalSummarizer, IntervalSummary};
// use crate::rada::modules::editfiltering::{ByMismatches, BySequencingQuality, Filters, FiltersBuilder};
use crate::rada::modules::workload::Workload;
use crate::rada::regions::ThreadContext;

// BySequencingQuality {mapq: 2, phread: 20};
// use rada::modules::stranding::{ByAtoIEditing, ByFeatures};
// use rada::per_region::ThreadContext;

// use crate::rada::other::LibraryType::Unstranded;
// use crate::rada::per_region::IntervalSummary;

mod rada;
// mod nuccounter;

// use rada::mismatches;
//
fn runall<Counter: NucCounter, Summarizer: IntervalSummarizer, Filter: IntervalSummaryFilter>(
    workload: Vec<Workload>,
    counter: Counter,
    summarizer: Summarizer,
    filter: Filter,
    bam: &Path,
) -> Vec<IntervalSummary> {
    let mut ctx = ThreadContext::new(&vec![bam], filter, counter, summarizer);

    workload
        .into_iter()
        .map(|w| rada::regions::run(&mut ctx, w))
        .flatten()
        .collect()

    // load bam file index only once in each thread
    // thread_local! {
    //     static CONTEXT: RefCell<Option<ThreadContext<Filter, Counter, Summarizer>>> = RefCell::new(None)
    // }
    // workload.into_par_iter()
    //     .map(move |w| {
    //         CONTEXT.with(|mut ctx| {
    //             if ctx.borrow().is_none() {
    //                 ctx.replace(Some(ThreadContext::new(&vec![bam], filter, counter, summarizer)));
    //             }
    //             let mut summaries = rada::regions::run(ctx.borrow_mut().as_mut().unwrap(), w);
    //             summaries
    //         })
    //     })
    //     .collect()
}

// fn runall(fasta: &Path,
//           workload: Vec<Workload>,
//           thresholds: &FilteringThresholds) {
//     // load bam file index only once in each thread
//     thread_local! {
//         static CONTEXT: RefCell<Option<Context>> = None
//     }
//     let bam = Path::new("/home/aleksander/projects/rada/tests/data/ADAR1-WT-IFNb-72h-Z22.bam");
//     for w in workload.into_iter() {
//         rayon::spawn(move || {
//             CONTEXT.with(|ctx| {
//                 if ctx.borrow().is_none() {
//                     ctx.replace(Some(Context::new(bam, fasta, 10_000)));
//                 }
//                 let mut ctx = ctx.borrow_mut().as_mut().unwrap();
//                 nuccounter::run(ctx, w, &thresholds);
//             })
//         })
//     }

// CONTEXT.set(&mut nuccounter::Context::new(bam, fasta, 10_000))
//
// workload.into_par_iter()
//     .map_init(
//         || CONTEXT.set(,
//         move |ctx, work| {
//
//         }
// );

// pool.scope(move |s| {
//     for w in workload {
//         let tx = tx.clone();
//         s.spawn(move |_| {
//             context.with(|ctx| {
//                 tx.send(nuccounter::run(ctx, w, thresholds)).expect("Failed to pass data between threads");
//             })
//         })
//     }
//     let mut a = 0;
//     for &_ in rx.into_iter() {
//         a += 1;
//     }
//     println!("{}", a);
// });
// }

fn main() {
    let settings = vec![
        ArgSettings::Required,
        ArgSettings::TakesValue,
        ArgSettings::AllowHyphenValues,
    ];

    let matches = App::new("rada")
        .args(vec![
            Arg::new("input")
                .short('i')
                .settings(&settings)
                .about("Bam file"),
            Arg::new("genome")
                .short('g')
                .settings(&settings)
                .about("Fasta file"),
            Arg::new("regions")
                .short('r')
                .settings(&settings)
                .about("Regions file"),
            Arg::new("annotation")
                .short('a')
                .settings(&settings)
                .about("Annotation gff3 file"),
            Arg::new("output")
                .short('o')
                .settings(&settings)
                .about("Output path"),
        ])
        .get_matches();

    let bam = PathBuf::from(matches.value_of("input").unwrap());
    let fasta = PathBuf::from(matches.value_of("genome").unwrap());
    assert!(bam.is_file() && fasta.is_file());

    // ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build_global()
    //     .expect("Failed to initialize global thread pool");
    // 1 - records regions
    let bed = PathBuf::from(matches.value_of("regions").unwrap());

    let output = matches.value_of("output").unwrap();
    let mut output = BufWriter::new(File::create(output).unwrap());

    writeln!(output, "chr\tstart\tend\tstrand\tname\t#A\t#T\t#G\t#C\tAA\tAG\tAC\tAT\tTA\tTG\tTC\tTT\tGA\tGG\tGC\tGT\tCA\tCG\tCC\tCT");

    let workload = Workload::from_bed_intervals(&bed, &bam);

    let gff3 = PathBuf::from(matches.value_of("annotation").unwrap());
    let strander: Vec<Box<dyn IntervalStrandPredictor>> = vec![
        Box::new(StrandByGenomicFeatures::from_gff3(&gff3)),
        Box::new(StrandByAtoIEditing::new(10, 0.05)),
    ];
    let strander = DynamicSequentialStrandPredictor::new(strander);

    let records_filter = RecordsFilterByQuality::new(2, 20);
    let maxsize = workload.iter().max_by_key(|x| x.len()).unwrap().len() as u32;
    let buffer = UnstrandedCountsBuffer::new(maxsize);
    let counter = BaseNucCounter::new(records_filter, buffer, Interval::new("".to_string(), 0..0));

    let refnuc = RefNucPredByHeurisitc::new(10, 0.95);

    let summary_filter = SummaryFilterByEditing::new(10, 0.05);

    let fasta = faidx::Reader::from_path(fasta).expect("Failed to parse Fasta");

    let summarizer = BaseSummarizer::new(fasta, refnuc, strander);

    for x in runall(workload, counter, summarizer, summary_filter, &bam) {
        // write!(output, "{}\t{}\t{}\t{}\t{}",
        //        x.interval.contig(), x.interval.range().start, x.interval.range().end, x.strand.strand_symbol(), x.name);
        // write!(output, "\t{}\t{}\t{}\t{}", x.refnuc.A, x.refnuc.T, x.refnuc.G, x.refnuc.C);
        //
        // for x in vec![&x.A, &x.T, &x.G, &x.C] {
        //     write!(output, "\t{}\t{}\t{}\t{}", x.A, x.G, x.C, x.T);
        // }
        write!(output, "\n");
    }
}
