use derive_more::Constructor;
use rust_htslib::bam::Record;

use super::{AlignedRead, ReadsFilter};

#[derive(Constructor)]
pub struct ReadsFilterByQuality {
    mapq: u8,
    phread: u8,
}

impl<R: AlignedRead> ReadsFilter<R> for ReadsFilterByQuality {
    fn is_read_ok(&self, record: &R) -> bool {
        record.mapq() != 255 && record.mapq() >= self.mapq
    }

    fn is_base_ok(&self, record: &R, base: usize) -> bool {
        record.qual()[base] >= self.phread
    }
    // 255 means no data is available
}
