use crate::core::dna::Nucleotide;

use bio_types::genome::Interval;
#[cfg(test)]
use mockall::{automock, predicate::*};

mod autoref;

use crate::core::dna::NucCounts;
pub use autoref::AutoRef;

#[cfg_attr(test, automock)]
pub trait RefEngine {
    fn reset(&mut self);
    fn run(&mut self, interval: &Interval, sequenced: &[NucCounts]);
    fn results<'a>(&'a self) -> Vec<&'a [Nucleotide]>;
}
