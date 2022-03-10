use std::convert::TryInto;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::{ReqStrand, Strand};
use flate2::read::GzDecoder;
use itertools::Itertools;

use crate::core::io;
use crate::core::mismatches::roi::{REATROIMismatchesVec, ROIMismatchesVec};
use crate::core::mismatches::site::{REATSiteMismatchesVec, SiteMismatchesVec};
use crate::core::mismatches::StrandingCounts;
use crate::core::stranding::predict::{StrandingAlgo, StrandingAlgoResult};

use super::utils;

#[derive(Clone)]
pub struct StrandByGenomicAnnotation {
    exons: AnnotMap<String, ReqStrand>,
    genes: AnnotMap<String, ReqStrand>,
}

impl StrandByGenomicAnnotation {
    pub fn from_gff(gff3: &Path, hook: impl Fn(&Contig<String, Strand>)) -> Self {
        io::utils::read_compressed!(gff3, Self::parse_gff, hook)
    }

    fn parse_gff<T: BufRead>(mut reader: T, hook: impl Fn(&Contig<String, Strand>)) -> Self {
        let mut exons: AnnotMap<String, ReqStrand> = AnnotMap::new();
        let mut genes: AnnotMap<String, ReqStrand> = AnnotMap::new();

        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to read annotation file") != 0 {
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
                    buf.clear();
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

    fn strand_in_index(&self, dummy: &Contig<String, Strand>, index: &AnnotMap<String, ReqStrand>) -> Strand {
        let intersected = index.find(dummy).unique_by(|x| x.data()).collect_vec();

        if intersected.len() == 1 {
            (*intersected.first().unwrap().data()).into()
        } else {
            Strand::Unknown
        }
    }

    fn predict(&self, contig: &str, range: Range<Position>) -> Strand {
        let (start, end) = (range.start, range.end);
        let dummy = Contig::new(contig.into(), start as isize, (end - start) as usize, Strand::Unknown);

        for index in [&self.exons, &self.genes] {
            let strand = self.strand_in_index(&dummy, index);
            if !strand.is_unknown() {
                return strand;
            }
        }
        Strand::Unknown
    }

    fn features_in(&self, supinterval: &impl AbstractInterval) -> Vec<Range<Position>> {
        let contig = supinterval.contig().to_owned();
        let (start, end) = (supinterval.range().start, supinterval.range().end);

        let key = Contig::new(contig, start as isize, (end - start) as usize, Strand::Unknown);
        let borders = (self.exons.find(&key))
            .chain(self.genes.find(&key))
            .map(|x| [x.interval().start as Position, x.interval().end as Position])
            .flatten()
            .sorted()
            .skip_while(|x| x <= &start)
            .dedup()
            .collect_vec();

        let mut supfeatures = Vec::with_capacity(borders.len() / 2 + 2);
        let mut prevind = start;
        for border in borders {
            let border = std::cmp::min(border, end);
            debug_assert!(border >= start);

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

impl StrandingAlgo<REATROIMismatchesVec> for StrandByGenomicAnnotation {
    fn predict(&self, x: &REATROIMismatchesVec) -> StrandingAlgoResult {
        let iter = x.postmasked().iter().map(|range| self.predict(x.contig(), range.clone()));
        utils::stranditer(iter)
    }
}

impl StrandingAlgo<REATSiteMismatchesVec> for StrandByGenomicAnnotation {
    fn predict(&self, x: &REATSiteMismatchesVec) -> StrandingAlgoResult {
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
        let mut counts = StrandingCounts::default();

        let mut iter = features.into_iter();
        let mut feature = iter.next().unwrap();
        let mut strand = None;
        for &pos in x.pos() {
            // While site is not inside the feature
            while !feature.contains(&pos) {
                feature = iter.next().unwrap();
                strand = None;
            }
            // Predict strand for the current feature (if not predicted yet)
            let strand = strand.unwrap_or_else(|| self.predict(x.contig(), feature.clone()));
            counts[strand] += 1;
            strands.push(strand);
        }

        StrandingAlgoResult::EachElement((strands, counts))
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use bio_types::strand::Same;

    use super::*;

    #[test]
    fn predict() {
        // create a dummy gff3 for the following genome (| - exons)
        // chr1: ---||||-->-|||||-->-||||-----
        // chr2: ---||-->---||-->--||-->------
        //       |||||||||----<-----|||--<----
        let gff3 = "\n\
        # I am comment number 1\n\
        # I am comment number 2\n\
        chr1\t.\tgene\t1\t29\t.\t+\t0\n\
        chr1\t.\texon\t4\t7\t.\t+\t0\n\
        chr1\t.\texon\t12\t16\t.\t+\t0\n\
        chr1\t.\texon\t21\t24\t.\t+\t0\n\
        # Forward gene\n\
        2\t.\tgene\t1\t29\t.\t+\t0\n\
        2\t.\texon\t4\t5\t.\t+\t0\n\
        2\t.\texon\t12\t13\t.\t+\t0\n\
        2\t.\texon\t19\t20\t.\t+\t0\n\
        # Reverse gene\n\
        2\t.\tgene\t1\t29\t.\t-\t0\n\
        2\t.\texon\t1\t9\t.\t-\t0\n\
        2\t.\texon\t20\t22\t.\t-\t0";

        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff3.as_bytes()), |_| {});
        for (contig, range, strand) in [
            ("chr1", 5..25, Strand::Forward),
            ("chr1", 3..4, Strand::Forward),
            ("chr1", 100..200, Strand::Unknown),
            ("chr1", 10..30, Strand::Forward),
            ("2", 0..3, Strand::Reverse),
            ("2", 3..8, Strand::Unknown),
            ("2", 9..10, Strand::Unknown),
            ("2", 12..13, Strand::Forward),
            ("2", 14..19, Strand::Forward),
            ("2", 21..50, Strand::Reverse),
            ("3", 1..5, Strand::Unknown),
        ] {
            let inferred = dummy.predict(contig, range);
            assert!(inferred.same(&strand))
        }
    }

    #[test]
    fn intervals_in() {
        let gff = "\n\
        chr1\t.\tgene\t2\t12\t.\t+\t0\n\
        chr1\t.\texon\t4\t6\t.\t+\t0\n\
        chr1\t.\texon\t9\t11\t.\t+\t0";
        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), |_| {});
        for (query, expected) in [
            (1..12, [1..3, 3..6, 6..8, 8..11, 11..12].to_vec()),
            (6..26, [6..8, 8..11, 11..12, 12..26].to_vec()),
            (2..7, [2..3, 3..6, 6..7].to_vec()),
            (2..7, [2..3, 3..6, 6..7].to_vec()),
            (0..7, [0..1, 1..3, 3..6, 6..7].to_vec()),
        ] {
            let inferred = dummy.features_in(&Interval::new("chr1".into(), query));
            assert_eq!(inferred, expected);
        }

        let gff = "
        chr1\t.\tgene\t11\t20\t.\t+\t0\n\
        chr1\t.\texon\t11\t12\t.\t+\t0\n\
        chr1\t.\texon\t17\t20\t.\t+\t0\n\
        #\n\
        chr1\t.\tgene\t17\t30\t.\t-\t0\n\
        chr1\t.\texon\t17\t24\t.\t-\t0\n\
        chr1\t.\texon\t29\t30\t.\t-\t0\n\
        #\n\
        chr1\t.\tgene\t2\t30\t.\t+\t0";
        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), |_| {});
        for (query, expected) in [
            (0..14, [0..1, 1..10, 10..12, 12..14].to_vec()),
            (13..30, [13..16, 16..20, 20..24, 24..28, 28..30].to_vec()),
            (16..36, [16..20, 20..24, 24..28, 28..30, 30..36].to_vec()),
        ] {
            let inferred = dummy.features_in(&Interval::new("chr1".into(), query));
            assert_eq!(inferred, expected);
        }
    }

    #[test]
    fn batched_sites() {
        // TODO: test for batched sites
    }
}
