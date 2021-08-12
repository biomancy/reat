use bio_types::genome::Interval;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use basecounter::BaseNucCounter;

use crate::rada::counting::CountsBufferContent;
use crate::rada::read::AlignedRead;

mod basecounter;

pub struct NucCounterContent<'a> {
    pub interval: Interval,
    pub counts: CountsBufferContent<'a>,
}

#[cfg_attr(test, automock)]
pub trait NucCounter<R: AlignedRead> {
    fn roi(&self) -> &Interval;
    fn process(&mut self, read: &mut R);
    fn reset(&mut self, roi: Interval);
    fn content(&'_ self) -> NucCounterContent<'_>;
}
