use std::ops::Range;
use std::path::PathBuf;

use bio_types::genome::Position;
use dyn_clone::DynClone;
#[cfg(test)]
use mockall::mock;
use rust_htslib::faidx;

use crate::core::dna::Nucleotide;

pub trait FastaReader: Send + DynClone {
    fn fetch(&mut self, contig: &str, range: Range<Position>);
    fn result(&self) -> &[Nucleotide];
}
dyn_clone::clone_trait_object!(FastaReader);

#[cfg(test)]
mock! {
    pub FastaReader {}
    impl Clone for FastaReader {
        fn clone(&self) -> Self;
    }

    impl FastaReader for FastaReader {
        fn fetch(&mut self, contig: &str, range: Range<Position>);
        fn result(&self) -> &[Nucleotide];
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
            faidx: faidx::Reader::from_path(&path).expect("Failed to open reference fasta file"),
            cache: Vec::new(),
            path,
        }
    }
}

impl FastaReader for BasicFastaReader {
    fn fetch(&mut self, contig: &str, range: Range<Position>) {
        self.cache.clear();

        let iter = self
            .faidx
            .fetch_seq(contig, range.start as usize, range.end as usize)
            .unwrap_or_else(|_| panic!("Failed to fetch sequence for region {}:{}-{}", contig, range.start, range.end))
            .iter()
            .map(|x| Nucleotide::from(*x));
        self.cache.extend(iter);

        let expected = (range.end - range.start) as usize;
        if self.cache.len() != expected {
            assert!(self.cache.len() > expected);
            self.cache.truncate(expected);
        }
        debug_assert_eq!(self.cache.len(), expected);
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
