
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
use crate::core::mismatches::roi::{ROIDataRef, ROIMismatchesVec};
use crate::core::mismatches::site::{SiteDataVec, SiteMismatchesVec};
use crate::core::mismatches::{MismatchesVec};
use crate::core::stranding::predict::StrandingAlgo;
use crate::core::strandutil::Stranded;

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

    fn features_in(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>> {
        let (start, end) = (range.start, range.end);

        let key = Contig::new(contig.to_owned(), start as isize, (end - start) as usize, Strand::Unknown);
        let borders = (self.exons.find(&key))
            .chain(self.genes.find(&key))
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
        let dummy = StrandByGenomicAnnotation::parse_gff(BufReader::new(gff.as_bytes()), |_| {});
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
    #[ignore]
    fn batched_sites() {
        todo!("Write tests for batch sites")
    }
}
