use std::ops::Range;

use bio_types::genome::Position;

pub use autoref::AutoRef;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::io::fasta::FastaReader;

mod autoref;
mod fasta;

pub struct RefEngineResult<'a> {
    pub predicted: &'a [Nucleotide],
    pub reference: &'a [Nucleotide],
}

pub trait RefEngine {
    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]);
    fn results(&self) -> RefEngineResult<'_>;
}
