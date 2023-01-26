use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::genome::Position;
use bio_types::strand::{ReqStrand, Strand};
use flate2::read::GzDecoder;
use itertools::Itertools;

use crate::core::io;
use crate::core::mismatches::roi::{ROIDataRef, ROIMismatchesVec};
use crate::core::mismatches::site::{SiteDataVec, SiteMismatchesVec};
use crate::core::mismatches::MismatchesVec;
use crate::core::stranding::predict::StrandingAlgo;
use crate::core::strandutil::Stranded;

use super::utils;

#[derive(Clone)]
pub struct StrandByGenomicAnnotation {
    exons: AnnotMap<String, ReqStrand>,
    genes: AnnotMap<String, ReqStrand>,
    extended3utr: AnnotMap<String, ReqStrand>,
}

impl StrandByGenomicAnnotation {
    pub fn from_gff(gff3: &Path, extended_3utr_size: u32, hook: impl Fn(usize)) -> Self {
        io::utils::read_compressed!(gff3, Self::parse_gff, extended_3utr_size, hook)
    }

    fn parse_gff<T: BufRead>(mut reader: T, extended_3utr_size: u32, hook: impl Fn(usize)) -> Self {
        let extended_3utr_size = extended_3utr_size as isize;

        let mut exons: AnnotMap<String, ReqStrand> = AnnotMap::new();
        let mut genes: AnnotMap<String, ReqStrand> = AnnotMap::new();
        let mut extended3utr: AnnotMap<String, ReqStrand> = AnnotMap::new();

        let mut parsedcnt: usize = 0;
        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to read annotation file") != 0 {
            if buf.starts_with('#') || buf == "\n" {
                buf.clear();
                continue;
            }

            let split: Vec<&str> = buf.split('\t').take(7).collect();
            assert_eq!(split.len(), 7);

            let strand = match split[6] {
                "+" => ReqStrand::Forward,
                "-" => ReqStrand::Reverse,
                _ => {
                    buf.clear();
                    continue;
                }
            };
            let (start, end) = (split[3].parse::<isize>().unwrap() - 1, split[4].parse::<isize>().unwrap());
            let record = Contig::new(split[0].into(), start, (end - start) as usize, Strand::Unknown);

            match split[2] {
                "exon" | "Exon" => exons.insert_at(strand, &record),
                "gene" | "Gene" => {
                    // 1 - insert into genes
                    genes.insert_at(strand, &record);
                    // 2 - infer extended utr size & insert into extended3utr
                    if extended_3utr_size > 0 {
                        let start = match strand {
                            ReqStrand::Forward => record.start() + record.length() as isize,
                            ReqStrand::Reverse => record.start() - extended_3utr_size,
                        };
                        let record = Contig::new(split[0].into(), start, extended_3utr_size as usize, Strand::Unknown);
                        extended3utr.insert_at(strand, &record);
                    }
                }
                _ => {
                    buf.clear();
                    continue;
                }
            };

            buf.clear();
            parsedcnt += 1;
            hook(parsedcnt);
        }
        StrandByGenomicAnnotation { exons, genes, extended3utr }
    }

    fn strand_in_index(&self, dummy: &Contig<String, Strand>, index: &AnnotMap<String, ReqStrand>) -> (u32, u32) {
        let (mut forward, mut reverse) = (0, 0);

        for req in index.find(dummy) {
            match req.data() {
                ReqStrand::Forward => forward += 1,
                ReqStrand::Reverse => reverse += 1,
            }
        }
        (forward, reverse)
    }

    fn predict(&self, contig: &str, range: Range<Position>) -> Strand {
        let (start, end) = (range.start, range.end);
        let dummy = Contig::new(contig.into(), start as isize, (end - start) as usize, Strand::Unknown);

        for index in [&self.exons, &self.genes, &self.extended3utr] {
            let (forward, reverse) = self.strand_in_index(&dummy, index);

            match (forward == 0, reverse == 0) {
                (true, true) => continue,                 // Nothing on both strands
                (true, false) => return Strand::Reverse,  // Features only on the forward strand
                (false, true) => return Strand::Forward,  // Features only on the reverse strand
                (false, false) => return Strand::Unknown, // Features on both strands
            };
        }
        Strand::Unknown
    }

    fn features_in(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>> {
        let (start, end) = (range.start, range.end);

        let key = Contig::new(contig.to_owned(), start as isize, (end - start) as usize, Strand::Unknown);
        let borders = (self.exons.find(&key))
            .chain(self.genes.find(&key))
            .chain(self.extended3utr.find(&key))
            .flat_map(|x| [x.interval().start as Position, x.interval().end as Position])
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

impl StrandingAlgo<ROIMismatchesVec> for StrandByGenomicAnnotation {
    fn predict(&self, contig: &str, items: &mut Stranded<ROIMismatchesVec>) {
        utils::assort_strands!(items, |x: ROIDataRef| self.predict(contig, x.roi.postmasked.clone()));
    }
}

impl StrandingAlgo<SiteMismatchesVec> for StrandByGenomicAnnotation {
    fn predict(&self, contig: &str, items: &mut Stranded<SiteMismatchesVec>) {
        if items.unknown.is_empty() {
            return;
        }
        let data = &mut items.unknown.data;
        let argsort = (0..data.len()).into_iter().sorted_by_key(|&x| data.pos[x]).collect_vec();

        // Get all annotated features in the given region (=regions with constant annotation)
        let range: Range<Position> = data.pos[*argsort.first().unwrap()]..data.pos[*argsort.last().unwrap()] + 1;
        let features = self.features_in(contig, range.clone());
        debug_assert!(
            !features.is_empty()
                && features.first().unwrap().start == range.start
                && features.last().unwrap().end == range.end
        );

        // Special case -> simply append all items to an existing vector
        if features.len() == 1 {
            match self.predict(contig, range) {
                Strand::Forward => items.forward.data.append(data),
                Strand::Reverse => items.reverse.data.append(data),
                Strand::Unknown => {}
            }
            return;
        }

        let mut remained = SiteDataVec::with_capacity(data.len() / 10);

        let mut iter = features.into_iter();
        let mut feature = iter.next().unwrap();
        let mut strand = self.predict(contig, feature.clone());
        for ind in argsort {
            let pos = data.pos[ind];
            // While site is not inside the feature
            while !feature.contains(&pos) {
                feature = iter.next().unwrap();

                // Predict strand if the next feature is useful
                if feature.contains(&pos) {
                    strand = self.predict(contig, feature.clone());
                }
            }

            let item = data.get(ind).unwrap();
            match strand {
                Strand::Forward => {
                    items.forward.data.push(item.into());
                }
                Strand::Reverse => {
                    items.reverse.data.push(item.into());
                }
                Strand::Unknown => remained.push(item.into()),
            }
        }
        items.unknown.data = remained;
    }
}

#[cfg(test)]
mod tests {
    use bio_types::strand::Same;
    use rand::seq::SliceRandom;
    use rand::thread_rng;

    use crate::core::mismatches::site::SiteData;

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

        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff3.as_bytes()), 0, |_| {});
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
            let inferred = dummy.predict(contig, range.clone());
            assert!(
                inferred.same(&strand),
                "{}:{}-{} ({}) vs ({})",
                contig,
                range.start,
                range.end,
                inferred.strand_symbol(),
                strand.strand_symbol()
            );
        }
    }

    #[test]
    fn intervals_in() {
        let gff = "\n\
        chr1\t.\tgene\t2\t12\t.\t+\t0\n\
        chr1\t.\texon\t4\t6\t.\t+\t0\n\
        chr1\t.\texon\t9\t11\t.\t+\t0";
        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), 0, |_| {});
        for (query, expected) in [
            (1..12, [1..3, 3..6, 6..8, 8..11, 11..12].to_vec()),
            (6..26, [6..8, 8..11, 11..12, 12..26].to_vec()),
            (2..7, [2..3, 3..6, 6..7].to_vec()),
            (2..7, [2..3, 3..6, 6..7].to_vec()),
            (0..7, [0..1, 1..3, 3..6, 6..7].to_vec()),
        ] {
            let inferred = dummy.features_in("chr1", query);
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
        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), 0, |_| {});
        for (query, expected) in [
            (0..14, [0..1, 1..10, 10..12, 12..14].to_vec()),
            (13..30, [13..16, 16..20, 20..24, 24..28, 28..30].to_vec()),
            (16..36, [16..20, 20..24, 24..28, 28..30, 30..36].to_vec()),
        ] {
            let inferred = dummy.features_in("chr1", query);
            assert_eq!(inferred, expected);
        }
    }

    #[test]
    fn batched_sites() {
        let gff = "\
        # Forward gene\n\
        1\tMyGene1\tgene\t4\t17\t0.132\t+\t0\n\
        1\tMyGene1\texon\t4\t5\t.\t+\t0\n\
        1\tMyGene1\texon\t9\t11\t1\t+\t0\n\
        1\tMyGene1\texon\t15\t17\t2\t+\t0\n\
        # Reverse gene\n\
        1\tMyGene2\tgene\t4\t18\t0\t-\t0\n\
        1\tMyGene2\texon\t4\t4\t.\t-\t0\n\
        1\tMyGene2\texon\t14\t16\t.\t-\t0\n\
        1\tMyGene2\texon\t18\t18\t.\t-\t0\n\
        ";

        let mut data = SiteDataVec::new();

        let mut range = (0..22).collect_vec();
        range.shuffle(&mut thread_rng());
        for pos in range {
            let mut d = SiteData::default();
            d.pos = pos;
            data.push(d);
        }
        let mm = SiteMismatchesVec::new("1".into(), Strand::Unknown, data);
        let mut workload = Stranded::with_fn(|strand| SiteMismatchesVec::new("1".into(), strand, SiteDataVec::new()));
        workload.unknown = mm;

        // Case 1 -> no extension
        let predictor = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), 0, |_| {});
        let mut m = workload.clone();
        StrandingAlgo::<SiteMismatchesVec>::predict(&predictor, "1", &mut m);

        debug_assert_eq!(m.forward.data.pos, &[4, 8, 9, 10, 16]);
        debug_assert_eq!(m.reverse.data.pos, &[13, 17]);
        debug_assert_eq!(m.unknown.data.pos, &[0, 1, 2, 3, 5, 6, 7, 11, 12, 14, 15, 18, 19, 20, 21]);

        // Case 2 -> has extension
        let predictor = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), 2, |_| {});
        let mut m = workload.clone();
        StrandingAlgo::<SiteMismatchesVec>::predict(&predictor, "1", &mut m);

        debug_assert_eq!(m.forward.data.pos, &[4, 8, 9, 10, 16, 18]);
        debug_assert_eq!(m.reverse.data.pos, &[1, 2, 13, 17]);
        debug_assert_eq!(m.unknown.data.pos, &[0, 3, 5, 6, 7, 11, 12, 14, 15, 19, 20, 21]);

        // Case 3 -> complete extension
        let predictor = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), 100, |_| {});
        let mut m = workload.clone();
        StrandingAlgo::<SiteMismatchesVec>::predict(&predictor, "1", &mut m);

        debug_assert_eq!(m.forward.data.pos, &[4, 8, 9, 10, 16, 18, 19, 20, 21]);
        debug_assert_eq!(m.reverse.data.pos, &[0, 1, 2, 13, 17]);
        debug_assert_eq!(m.unknown.data.pos, &[3, 5, 6, 7, 11, 12, 14, 15]);
    }
}
