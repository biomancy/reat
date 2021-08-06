// use std::path::{Path, PathBuf};
// use std::collections::HashMap;
//
// use bio_types::strand::ReqStrand;
// use fstrings::{f, format_args_f};
// use rust_htslib::{bam, bam::Read, faidx};
//
// use super::modules::{filtering, counting};
// use super::modules::counting::{Counts, CountsBuffer};
// use super::modules::workload::Workload;
// use std::fs::records;
// use crate::rada::other::LibraryType;
// use bio_types::genome::AbstractInterval;
//
// #[allow(non_snake_case)]
// pub struct MismatchesSummary {
//     // TODO: move to the outer level?
//     pub tid: u32,
//     pub pos: u32,
//     pub plscounts: Counts,
//     pub mnscounts: Counts,
//     pub reference: u8,
// }
//
// pub struct ThreadContext {
//     refreader: faidx::Reader,
//     htsreaders: HashMap<PathBuf, bam::IndexedReader>,
//     buffer: CountsBuffer
// }
//
//
// impl ThreadContext {
//     pub fn new(htsfiles: &[&Path], reference: &Path, max_worksize: u32) -> Self {
//         let mut htsreaders: HashMap<PathBuf, bam::IndexedReader> = HashMap::new();
//         for &hts in htsfiles {
//             let reader = bam::IndexedReader::from_path(&hts)
//                 .expect(&f!("Failed to open file {file}", file=hts.display()));
//             htsreaders.insert(PathBuf::from(hts), reader);
//         }
//         let refreader = faidx::Reader::from_path(reference)
//             .expect(&f!("Failed to open file {file}", file=reference.display()));
//
//         let buffer = CountsBuffer::new(max_worksize);
//         println!("Created");
//         ThreadContext { refreader, htsreaders, buffer }
//     }
// }
//
// #[inline]
// pub fn reverse_complement(nuc: &u8) -> u8 {
//     match nuc {
//         b'A' | b'a' => b'T',
//         b'T' | b't' => b'A',
//         b'G' | b'g' => b'C',
//         b'C' | b'c' => b'G',
//         _ => *nuc
//     }
// }
//
//
// pub fn run(ctx: &mut ThreadContext, library: LibraryType,
//            work: Workload, filters: &filtering::Filters) -> Vec<MismatchesSummary> {
//     // 1. reset the counting
//     ctx.buffer.reset();
//
//     // 2. do counting
//     let mut htsreader = ctx.htsreaders.get_mut(work.bam)
//         .expect("Bug: thread failed to initialize all BAM readers");
//     counting::in_region(htsreader, library, &work.interval, &mut ctx.buffer, filters.by_quality());
//
//     // 3. fetch reference sequence
//     let (start, end) = (work.interval.range().start as usize, work.interval.range().end as usize - 1);
//     let reference = ctx.refreader.fetch_seq(work.interval.contig(), start, end)
//         .expect(&f!("Failed to fetch sequence for region TODO:{start}-{end}"));
//
//     // 4. Auto-skip obvious SNPs
//     let worksize = work.interval.range().end as usize - work.interval.range().start as usize;
//     let mut summarization = Vec::with_capacity(worksize / 2);
//     debug_assert!(reference.len() == worksize);
//
//     let autorefthr = 1.0f32 - 0.001f32;
//     let minmis = filters.by_mismatches().as_ref().unwrap();
//     for (ind, mut refnuc) in reference.iter().enumerate() {
//         let plscounts = ctx.buffer.plstrand()[ind];
//         let mnscounts = ctx.buffer.mnstrand()[ind];
//
//         let coverage = plscounts.coverage() + mnscounts.coverage();
//         if coverage == 0 {
//             continue;
//         }
//         // TODO: abstract this
//         // TODO: I need to use here some fun algo:
//         // https://en.wikipedia.org/wiki/SNV_calling_from_NGS_data
//         if coverage > 10 {
//             if (plscounts.A + mnscounts.T) as f32 / coverage as f32  > autorefthr {
//                 refnuc = &b'A';
//             } else if (plscounts.T + mnscounts.A) as f32 / coverage as f32  > autorefthr {
//                 refnuc = &b'T';
//             } else if (plscounts.G + mnscounts.C) as f32 / coverage as f32  > autorefthr {
//                 refnuc = &b'G';
//             } else if (plscounts.C + mnscounts.G) as f32 / coverage as f32  > autorefthr {
//                 refnuc = &b'C';
//             };
//         }
//
//         let coverage = coverage as f32;
//         let mismatches = plscounts.mismatches(refnuc) + mnscounts.mismatches(refnuc);
//         let isok = (mismatches >= minmis.number) &&
//             (mismatches as f32 / coverage > minmis.freq);
//         if !isok {
//             continue;
//         }
//
//         // Infer strand here
//
//         summarization.push(
//             MismatchesSummary{
//                 tid: work.interval.tid,
//                 pos: (start + ind) as u32,
//                 plscounts,
//                 mnscounts,
//                 reference: *refnuc
//             }
//         )
//     }
//     summarization
// }
