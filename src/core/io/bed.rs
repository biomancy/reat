use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use bio_types::genome::{AbstractInterval, Interval, Position};
use derive_getters::Dissolve;
use flate2::bufread::GzDecoder;
use mockall::Any;

#[derive(Eq, PartialEq, Debug, Dissolve)]
pub struct BedRecord {
    pub name: String,
    pub interval: Interval,
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
        let split: Vec<&str> = line.split('\t').take(4).collect();
        assert!(split.len() >= 3);

        let start = split[1].parse().expect("Failed to filters string start");
        let end = split[2].parse().expect("Failed to filters string start");
        let interval = Interval::new(split[0].to_owned(), Range { start, end });

        let name = if split.len() >= 4 { split[3] } else { &"" };
        let name = name.to_string();

        records.push(BedRecord { name, interval });
        buf.clear();
    }
    records
}

pub fn parse(bed: impl AsRef<Path>) -> Vec<BedRecord> {
    let bed = bed.as_ref();
    let extensions = bed.file_name().and_then(OsStr::to_str).expect(
        "Failed to infer extension for BED file with target regions. \
                          Only \".bed\" and \".bed.gz\" files are supported",
    );
    let extensions: Vec<&str> = extensions.split('.').collect();
    assert!(!extensions.is_empty());

    let file = File::open(bed).expect("Failed to open bed file");
    let file = BufReader::new(file);
    let last = *extensions.last().unwrap();

    if last == "gz" {
        assert_eq!(extensions[extensions.len() - 2], "bed");
        _parse(BufReader::new(GzDecoder::new(file)))
    } else {
        assert_eq!(last, "bed");
        _parse(file)
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use bio_types::genome::Position;

    use super::*;

    fn br(chr: &str, range: Range<Position>, name: &str) -> BedRecord {
        BedRecord { interval: Interval::new(chr.to_string(), range), name: name.to_string() }
    }

    #[test]
    fn empty() {
        let bed = "";
        assert!(_parse(BufReader::new(bed.as_bytes())).is_empty());
    }

    #[test]
    fn correct() {
        let bed = "\
        chr1\t10\t20\tReg1\n\
        chr1\t50\t60\tIII\n\
        chr1\t30\t40\t2\n\
        chr1\t70\t80\t\n";
        let records = vec![
            br("chr1", 10..20, "Reg1"),
            br("chr1", 50..60, "III"),
            br("chr1", 30..40, "2"),
            br("chr1", 70..80, ""),
        ];

        assert_eq!(records, _parse(BufReader::new(bed.as_bytes())));
    }

    #[test]
    fn empty_lines() {
        let bed = "\
        Very-long-line\t1000000\t2000000\tLorem_Ipsum_Doler_Sit_A_Met\t\t\t\n
        \n\
        \n\
        MT\t10\t12\t1R1\n\
        \n\
        1\t30\t301\t.\n\
        chr4\t700\t1800\t\"1\"\n
        \n\
        \n";
        let records = vec![
            br("Very-long-line", 1000000..2000000, "Lorem_Ipsum_Doler_Sit_A_Met"),
            br("MT", 10..12, "1R1"),
            br("1", 30..301, "."),
            br("chr4", 700..1800, "\"1\""),
        ];

        assert_eq!(records, _parse(BufReader::new(bed.as_bytes())));
    }
}
