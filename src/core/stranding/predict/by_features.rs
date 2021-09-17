use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::genome::{AbstractInterval, AbstractLocus, Interval, Locus};
use bio_types::strand::{ReqStrand, Strand};
use flate2::bufread::GzDecoder;
use itertools::Itertools;

use crate::core::counting::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::summary::MismatchesSummary;

use super::{LocusStrandPredictor, ROIStrandPredictor, StrandPredictor};

#[derive(Clone)]
pub struct StrandByGenomicFeatures {
    exons: AnnotMap<String, ReqStrand>,
    genes: AnnotMap<String, ReqStrand>,
}

impl StrandPredictor for StrandByGenomicFeatures {}

impl StrandByGenomicFeatures {
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
            StrandByGenomicFeatures::parse_gff3(io::BufReader::new(GzDecoder::new(gff3)), hook)
        } else {
            assert_eq!(last, "gff3");
            StrandByGenomicFeatures::parse_gff3(gff3, hook)
        }
    }

    fn parse_gff3<T: BufRead>(mut reader: T, hook: impl Fn(&Contig<String, Strand>)) -> Self {
        let mut exons: AnnotMap<String, ReqStrand> = AnnotMap::new();
        let mut genes: AnnotMap<String, ReqStrand> = AnnotMap::new();

        let mut buf = String::new();
        while reader.read_line(&mut buf).expect("Failed to reads GFF3 file") != 0 {
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
        StrandByGenomicFeatures { exons, genes }
    }

    fn predict_by_index(&self, interval: &Interval, index: &AnnotMap<String, ReqStrand>) -> Strand {
        let (start, end) = (interval.range().start, interval.range().end);
        let dummy = Contig::new(interval.contig().into(), start as isize, (end - start) as usize, Strand::Unknown);
        let intersected = index.find(&dummy).unique_by(|x| x.data()).collect_vec();
        if intersected.len() == 1 {
            Strand::from(*intersected[0].data())
        } else {
            Strand::Unknown
        }
    }

    pub fn infer_strand(&self, interval: &Interval) -> Strand {
        let byexons = self.predict_by_index(interval, &self.exons);
        if byexons.is_unknown() {
            self.predict_by_index(interval, &self.genes)
        } else {
            byexons
        }
    }
}

impl ROIStrandPredictor for StrandByGenomicFeatures {
    fn predict(&self, interval: &Interval, _: &MismatchesSummary) -> Strand {
        self.infer_strand(interval)
    }
}

impl LocusStrandPredictor for StrandByGenomicFeatures {
    fn predict(&self, locus: &Locus, _: &Nucleotide, _: &NucCounts) -> Strand {
        self.infer_strand(&Interval::new(locus.contig().to_string(), locus.pos()..locus.pos() + 1))
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
        # Forward gene
        2\t.\tgene\t1\t29\t.\t+\t0\n\
        2\t.\texon\t4\t5\t.\t+\t0\n\
        2\t.\texon\t12\t13\t.\t+\t0\n\
        2\t.\texon\t19\t20\t.\t+\t0\n\
        # Reverse gene
        2\t.\tgene\t1\t29\t.\t-\t0\n\
        2\t.\texon\t1\t9\t.\t-\t0\n\
        2\t.\texon\t20\t22\t.\t-\t0";

        let dummy = StrandByGenomicFeatures::parse_gff3(io::BufReader::new(gff3.as_bytes()), |_| {});
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
            let inferred = dummy.infer_strand(&Interval::new(contig.into(), range));
            assert!(inferred.same(&strand))
        }
    }
}
