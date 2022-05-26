use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::{Same, Strand};
use derive_getters::Dissolve;
use flate2::bufread::GzDecoder;

use super::utils;

#[derive(Debug, Clone, Dissolve)]
pub struct BedRecord {
    pub name: String,
    pub strand: Strand,
    pub interval: Interval,
}

impl PartialEq for BedRecord {
    fn eq(&self, other: &Self) -> bool {
        self.strand.same(&other.strand) && self.name == other.name && self.interval == other.interval
    }
}

impl AbstractInterval for BedRecord {
    fn contig(&self) -> &str {
        self.interval.contig()
    }

    fn range(&self) -> Range<Position> {
        self.interval.range()
    }
}

fn _parse<T: BufRead>(mut reader: T) -> Vec<BedRecord> {
    let mut records = Vec::new();

    let mut buf = String::new();
    while reader.read_line(&mut buf).expect("Failed to read BED file") != 0 {
        let line = buf.trim_end();
        if line.is_empty() {
            buf.clear();
            continue;
        }
        let split: Vec<&str> = line.split('\t').take(6).collect();
        assert!(split.len() >= 3);

        let start = split[1].parse().expect("Failed to filters string start");
        let end = split[2].parse().expect("Failed to filters string start");
        assert!(end > start, "{}", line);
        let interval = Interval::new(split[0].to_owned(), Range { start, end });

        let name = split.get(3).unwrap_or(&"").to_string();
        let strand = split.get(5).map_or(Strand::Unknown, |x| {
            Strand::from_char(&x.chars().next().unwrap()).expect("Failed to parse strand")
        });

        records.push(BedRecord { name, strand, interval });
        buf.clear();
    }
    records
}

pub fn parse(bed: impl AsRef<Path>) -> Vec<BedRecord> {
    let bed = bed.as_ref();
    utils::read_compressed!(bed, _parse)
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use bio_types::genome::Position;

    use super::*;

    fn br(chr: &str, range: Range<Position>, name: &str, strand: Strand) -> BedRecord {
        BedRecord { interval: Interval::new(chr.to_string(), range), name: name.to_string(), strand }
    }

    #[test]
    fn empty() {
        let bed = "";
        assert!(_parse(BufReader::new(bed.as_bytes())).is_empty());
    }

    #[test]
    fn correct() {
        let bed = "\
        chr1\t10\t20\tReg1\t.\t.\n\
        chr1\t50\t60\tIII\t1\t.\n\
        chr1\t30\t40\t2\t.\t-\n\
        chr1\t70\t80\t\t.\t+\n";
        let records = vec![
            br("chr1", 10..20, "Reg1", Strand::Unknown),
            br("chr1", 50..60, "III", Strand::Unknown),
            br("chr1", 30..40, "2", Strand::Reverse),
            br("chr1", 70..80, "", Strand::Forward),
        ];

        assert_eq!(records, _parse(BufReader::new(bed.as_bytes())));
    }

    #[test]
    fn empty_lines() {
        let bed = "\
        Very-long-line\t1000000\t2000000\tLorem_Ipsum_Doler_Sit_A_Met\t\t+\t\n
        \n\
        \n\
        MT\t10\t12\t1R1\t113\t.\n\
        \n\
        1\t30\t301\t.\n\
        chr4\t700\t1800\t\"1\"\t3\t-\n
        \n\
        \n";
        let records = vec![
            br("Very-long-line", 1000000..2000000, "Lorem_Ipsum_Doler_Sit_A_Met", Strand::Forward),
            br("MT", 10..12, "1R1", Strand::Unknown),
            br("1", 30..301, ".", Strand::Unknown),
            br("chr4", 700..1800, "\"1\"", Strand::Reverse),
        ];

        assert_eq!(records, _parse(BufReader::new(bed.as_bytes())));
    }
}
