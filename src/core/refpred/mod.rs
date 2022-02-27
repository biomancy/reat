use std::ops::Range;

use bio_types::genome::{Interval, Position};
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use autoref::AutoRef;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

mod autoref;

#[cfg_attr(test, automock)]
pub trait RefEngine {
    fn reset(&mut self);
    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]);
    fn results<'a>(&'a self) -> &'a [Nucleotide];
}
