use std::path::PathBuf;

use bio_types::genome::{AbstractInterval, Interval};
use derive_more::Constructor;
#[cfg(test)]
use mockall::mock;
use rust_htslib::faidx;

use crate::core::dna::Nucleotide;

pub trait FastaReader {
    fn fetch(&mut self, interval: &Interval);
    fn result(&self) -> &[Nucleotide];
}

#[cfg(test)]
mock! {
    pub FastaReader {}
    impl FastaReader for FastaReader {
        fn fetch(&mut self, interval: &Interval);
        fn result(&self) -> &[Nucleotide];
    }
    impl Clone for FastaReader {
        fn clone(&self) -> Self;
    }
}

pub struct BasicFastaReader {
    faidx: faidx::Reader,
    cache: Vec<Nucleotide>,
    path: PathBuf,
}
unsafe impl Send for BasicFastaReader {}
unsafe impl Sync for BasicFastaReader {}

impl BasicFastaReader {
    pub fn new(path: PathBuf) -> Self {
        Self {
            faidx: faidx::Reader::from_path(&path).expect("Failed to open requested fasta file"),
            cache: Vec::new(),
            path,
        }
    }
}

impl FastaReader for BasicFastaReader {
    fn fetch(&mut self, interval: &Interval) {
        self.cache.clear();

        let iter = self
            .faidx
            .fetch_seq(interval.contig(), interval.range().start as usize, interval.range().end as usize)
            .unwrap_or_else(|_| panic!("Failed to fetch sequence for region {:?}", interval))
            .into_iter()
            .map(|x| Nucleotide::from(*x));
        self.cache.extend(iter);
    }

    fn result(&self) -> &[Nucleotide] {
        &self.cache
    }
}

impl Clone for BasicFastaReader {
    fn clone(&self) -> Self {
        Self::new(self.path.clone())
    }
}
