use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::ops::Range;
use std::path::Path;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::genome::{AbstractInterval, AbstractLocus, Interval, Position};
use bio_types::strand::{ReqStrand, Strand};
use flate2::bufread::GzDecoder;
use itertools::Itertools;

use crate::core::mismatches::roi::{BatchedROIMismatches, REATBatchedROIMismatches};
use crate::core::mismatches::site::{BinnedSiteMismatches, REATBatchedSiteMismatches};
use crate::core::mismatches::BatchedMismatches;
use crate::core::stranding::predict::ctx::StrandingAlgoResult;
use crate::core::stranding::predict::{StrandingAlgo, StrandingContext};

#[derive(Clone)]
pub struct StrandByGenomicAnnotation {
    exons: AnnotMap<String, ReqStrand>,
    genes: AnnotMap<String, ReqStrand>,
}

impl StrandByGenomicAnnotation {
    pub fn from_gff3(gff3: &Path, hook: impl Fn(&Contig<String, Strand>)) -> Self {
        let extensions = gff3
            .file_name()
            .and_then(OsStr::to_str)
            .expect("Failed to infer extension for regions file. Only .bed and .bed.gz are supported");
        let extensions: Vec<&str> = extensions.split('.').collect();
        assert!(!extensions.is_empty());

        let gff3 = File::open(gff3).expect("Failed to open GFF3 file");
        let gff3 = io::BufReader::new(gff3);
        let last = *extensions.last().unwrap();

        if last == "gz" || last == "gzip" {
            assert_eq!(extensions[extensions.len() - 2], "gff3");
            StrandByGenomicAnnotation::parse_gff3(io::BufReader::new(GzDecoder::new(gff3)), hook)
        } else {
            assert_eq!(last, "gff3");
            StrandByGenomicAnnotation::parse_gff3(gff3, hook)
        }
    }

    fn parse_gff3<T: BufRead>(mut reader: T, hook: impl Fn(&Contig<String, Strand>)) -> Self {
        let mut exons: AnnotMap<String, ReqStrand> = AnnotMap::new();
        let mut genes: AnnotMap<String, ReqStrand> = AnnotMap::new();

        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to filters GFF3 file") != 0 {
            if buf.starts_with('#') || buf == "\n" {
                buf.clear();
                continue;
            }

            let split: Vec<&str> = buf.split('\t').take(7).collect();
            assert_eq!(split.len(), 7);

            let interval_tree = match split[2] {
                "exon" => &mut exons,
                "gene" => &mut genes,
                _ => {
                    buf.clear();
                    continue;
                }
            };

            let (start, end) = (split[3].parse::<isize>().unwrap() - 1, split[4].parse::<isize>().unwrap());
            let strand = match split[6] {
                "+" => ReqStrand::Forward,
                "-" => ReqStrand::Reverse,
                _ => {
                    continue;
                }
            };
            let record = Contig::new(split[0].into(), start, (end - start) as usize, Strand::Unknown);
            hook(&record);
            interval_tree.insert_at(strand, &record);

            buf.clear();
        }
        StrandByGenomicAnnotation { exons, genes }
    }

    fn predict_by_index(&self, dummy: &Contig<String, Strand>, index: &AnnotMap<String, ReqStrand>) -> Strand {
        let intersected = index.find(dummy).unique_by(|x| x.data()).collect_vec();
        if intersected.len() == 1 {
            Strand::from(*intersected[0].data())
        } else {
            Strand::Unknown
        }
    }

    pub fn predict(&self, contig: &str, range: Range<Position>) -> Strand {
        let (start, end) = (range.start, range.end);
        let dummy = Contig::new(contig.into(), start as isize, (end - start) as usize, Strand::Unknown);

        let byexons = self.predict_by_index(&dummy, &self.exons);
        if byexons.is_unknown() {
            self.predict_by_index(&dummy, &self.genes)
        } else {
            byexons
        }
    }

    pub fn features_in(&self, supinterval: &impl AbstractInterval) -> Vec<Range<Position>> {
        let contig = supinterval.contig().to_owned();
        let (start, end) = (supinterval.range().start, supinterval.range().end);

        let key = Contig::new(contig, start as isize, (end - start) as usize, Strand::Unknown);
        let borders = (self.exons.find(&key))
            .chain(self.genes.find(&key))
            .map(|x| [x.interval().start as u64, x.interval().end as u64])
            .flatten()
            .sorted()
            .skip_while(|x| x <= &start)
            .dedup()
            .collect_vec();

        let mut supfeatures = Vec::with_capacity(borders.len() / 2 + 2);
        let mut prevind = start;
        for border in borders {
            let border = std::cmp::min(border, end);
            supfeatures.push(prevind..border);
            if border >= end {
                break;
            }
            prevind = border;
        }
        if supfeatures.is_empty() || supfeatures.last().unwrap().end < end {
            supfeatures.push(prevind..end);
        }
        debug_assert!(!supfeatures.is_empty());
        debug_assert!(supfeatures.windows(2).all(|x| x[0].end == x[1].start));
        debug_assert!(supfeatures.first().unwrap().start == start);
        debug_assert!(supfeatures.last().unwrap().end == end);
        supfeatures
    }
}

impl StrandingAlgo<REATBatchedROIMismatches> for StrandByGenomicAnnotation {
    fn predict(
        &self,
        mut ctx: StrandingContext<REATBatchedROIMismatches>,
    ) -> StrandingContext<REATBatchedROIMismatches> {
        ctx.apply(|x| {
            StrandingAlgoResult::EachElement(
                x.rois().iter().map(|range| self.predict(x.contig(), range.clone())).collect(),
            )
        });
        ctx
    }
}

impl StrandingAlgo<REATBatchedSiteMismatches> for StrandByGenomicAnnotation {
    fn predict(
        &self,
        mut ctx: StrandingContext<REATBatchedSiteMismatches>,
    ) -> StrandingContext<REATBatchedSiteMismatches> {
        ctx.apply(|x| {
            // Get all annotated features in the given region (=regions with constant annotation)
            let features = self.features_in(x);
            debug_assert!(
                features.len() >= 1
                    && features.first().unwrap().start == x.range().start
                    && features.last().unwrap().end == x.range().end
            );

            if features.len() == 1 {
                let strand = self.predict(x.contig(), x.range());
                return StrandingAlgoResult::AllElements(strand);
            }

            let mut strands = Vec::with_capacity(x.pos().len());

            let mut iter = features.into_iter();
            let mut feature = iter.next().unwrap();
            let mut strand = None;
            for &pos in x.pos() {
                // While site is not inside the feature
                while !(feature.start <= pos && pos < feature.end) {
                    feature = iter.next().unwrap();
                    strand = None;
                }
                // Predict strand for the current feature
                if strand.is_none() {
                    strand = Some(self.predict(x.contig(), feature.clone()));
                }
                strands.push(strand.unwrap());
            }

            StrandingAlgoResult::EachElement(strands)
        });
        ctx
    }
}

// impl<T: IntermediateROIMismatches> ROIStrandPredictor<T> for StrandByGenomicAnnotation {
//     fn predict(&self, mut rois: Vec<T>) -> Vec<T> {
//         for r in &mut rois {
//             r.set_strand(self.predict(r.roi()));
//         }
//         rois
//     }
// }
//
// impl<T: IntermediateIntervalMismatches> IntervalStrandPredictor<T> for StrandByGenomicAnnotation {
//     fn predict(&self, blocks: Vec<T>) -> Vec<T> {
//         if blocks.is_empty() {
//             return blocks;
//         }
//
//         let mut result = Vec::with_capacity(blocks.len());
//         for mut i in blocks {
//             debug_assert!(i.strand().is_unknown());
//
//             // Split region into sub intervals with constant annotation if needed
//             let interval = i.interval();
//             let supfeatures = self.intervals_in(interval);
//             debug_assert!(
//                 supfeatures.len() >= 1
//                     && supfeatures.first().unwrap().range().start == i.interval().range().start
//                     && supfeatures.last().unwrap().range().end == i.interval().range().end
//             );
//
//             if supfeatures.len() == 1 {
//                 let strand = self.predict(interval);
//                 i.set_strand(strand);
//                 result.push(i);
//                 continue;
//             }
//
//             let interst = interval.range().start;
//             let splits = supfeatures.iter().skip(1).map(|x| (x.range().start - interst) as usize).collect_vec();
//             let splited = i.split(&splits);
//
//             debug_assert!(splited.iter().map(|x| x.interval()).zip(supfeatures.iter()).all(|(x, y)| x == y));
//
//             for mut b in splited {
//                 b.set_strand(self.predict(b.interval()));
//                 result.push(b);
//             }
//         }
//         result
//     }
// }

// #[cfg(test)]
// mod tests {
//     use bio_types::strand::Same;
//
//     use super::*;
//
//     #[test]
//     fn predict() {
//         // create a dummy gff3 for the following genome (| - exons)
//         // chr1: ---||||-->-|||||-->-||||-----
//         // chr2: ---||-->---||-->--||-->------
//         //       |||||||||----<-----|||--<----
//         let gff3 = "\n\
//         # I am comment number 1\n\
//         # I am comment number 2\n\
//         chr1\t.\tgene\t1\t29\t.\t+\t0\n\
//         chr1\t.\texon\t4\t7\t.\t+\t0\n\
//         chr1\t.\texon\t12\t16\t.\t+\t0\n\
//         chr1\t.\texon\t21\t24\t.\t+\t0\n\
//         # Forward gene\n\
//         2\t.\tgene\t1\t29\t.\t+\t0\n\
//         2\t.\texon\t4\t5\t.\t+\t0\n\
//         2\t.\texon\t12\t13\t.\t+\t0\n\
//         2\t.\texon\t19\t20\t.\t+\t0\n\
//         # Reverse gene\n\
//         2\t.\tgene\t1\t29\t.\t-\t0\n\
//         2\t.\texon\t1\t9\t.\t-\t0\n\
//         2\t.\texon\t20\t22\t.\t-\t0";
//
//         let dummy = StrandByGenomicAnnotation::parse_gff3(io::BufReader::new(gff3.as_bytes()), |_| {});
//         for (contig, range, strand) in [
//             ("chr1", 5..25, Strand::Forward),
//             ("chr1", 3..4, Strand::Forward),
//             ("chr1", 100..200, Strand::Unknown),
//             ("chr1", 10..30, Strand::Forward),
//             ("2", 0..3, Strand::Reverse),
//             ("2", 3..8, Strand::Unknown),
//             ("2", 9..10, Strand::Unknown),
//             ("2", 12..13, Strand::Forward),
//             ("2", 14..19, Strand::Forward),
//             ("2", 21..50, Strand::Reverse),
//             ("3", 1..5, Strand::Unknown),
//         ] {
//             let inferred = dummy.predict(&Interval::new(contig.into(), range));
//             assert!(inferred.same(&strand))
//         }
//     }
//
//     #[test]
//     fn intervals_in() {
//         let interval = |range| Interval::new("chr1".into(), range);
//
//         let gff3 = "\n\
//         chr1\t.\tgene\t2\t12\t.\t+\t0\n\
//         chr1\t.\texon\t4\t6\t.\t+\t0\n\
//         chr1\t.\texon\t9\t11\t.\t+\t0";
//         let dummy = StrandByGenomicAnnotation::parse_gff3(io::BufReader::new(gff3.as_bytes()), |_| {});
//         for (query, expected) in [
//             (1..12, [1..3, 3..6, 6..8, 8..11, 11..12].to_vec()),
//             (6..26, [6..8, 8..11, 11..12, 12..26].to_vec()),
//             (2..7, [2..3, 3..6, 6..7].to_vec()),
//             (2..7, [2..3, 3..6, 6..7].to_vec()),
//             (0..7, [0..1, 1..3, 3..6, 6..7].to_vec()),
//         ] {
//             let expected: Vec<Interval> = expected.into_iter().map(|x| interval(x)).collect();
//             let inferred = dummy.intervals_in(&interval(query));
//             assert_eq!(inferred, expected);
//         }
//
//         let gff3 = "
//         chr1\t.\tgene\t11\t20\t.\t+\t0\n\
//         chr1\t.\texon\t11\t12\t.\t+\t0\n\
//         chr1\t.\texon\t17\t20\t.\t+\t0\n\
//         #\n\
//         chr1\t.\tgene\t17\t30\t.\t-\t0\n\
//         chr1\t.\texon\t17\t24\t.\t-\t0\n\
//         chr1\t.\texon\t29\t30\t.\t-\t0\n\
//         #\n\
//         chr1\t.\tgene\t2\t30\t.\t+\t0";
//         let dummy = StrandByGenomicAnnotation::parse_gff3(io::BufReader::new(gff3.as_bytes()), |_| {});
//         for (query, expected) in [
//             (0..14, [0..1, 1..10, 10..12, 12..14].to_vec()),
//             (13..30, [13..16, 16..20, 20..24, 24..28, 28..30].to_vec()),
//             (16..36, [16..20, 20..24, 24..28, 28..30, 30..36].to_vec()),
//         ] {
//             let expected: Vec<Interval> = expected.into_iter().map(|x| interval(x)).collect();
//             let inferred = dummy.intervals_in(&interval(query));
//             assert_eq!(inferred, expected);
//         }
//     }
// }
