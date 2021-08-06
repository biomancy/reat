use rust_htslib::bam::Record;

mod byquality;

pub use byquality::RecordsFilterByQuality;

pub trait RecordsFilter {
    fn is_record_ok(&self, record: &Record) -> bool;
    fn is_base_ok(&self, record: &Record, qual: &[u8], base: usize) -> bool;
}
