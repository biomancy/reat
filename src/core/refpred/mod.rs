use std::ops::Range;

use bio_types::genome::{Interval, Position};
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use autoref::AutoRef;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

mod autoref;

pub struct RefEngineResult<'a> {
    pub predicted: &'a [Nucleotide],
    pub reference: &'a [Nucleotide],
}

#[cfg_attr(test, automock)]
pub trait RefEngine {
    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]);
    fn results(&self) -> RefEngineResult<'_>;
}
