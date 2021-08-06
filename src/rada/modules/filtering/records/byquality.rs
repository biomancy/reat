use bio_types::sequence::SequenceRead;
use rust_htslib::bam::Record;

use super::RecordsFilter;

use derive_more::Constructor;

#[derive(Constructor)]
pub struct RecordsFilterByQuality {
    mapq: u8,
    phread: u8,
}

impl RecordsFilter for RecordsFilterByQuality {
    // 255 means no data is available
    fn is_record_ok(&self, record: &Record) -> bool {
        record.mapq() != 255 && record.mapq() >= self.mapq
    }

    fn is_base_ok(&self, _: &Record, qual: &[u8], base: usize) -> bool {
        qual[base] >= self.phread
    }
}
