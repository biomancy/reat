use bio_types::genome::Interval;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use basecounter::BaseNucCounter;

use crate::core::counting::CountsBufferContent;
use crate::core::read::AlignedRead;

mod basecounter;

#[derive(Clone)]
pub struct NucCounterContent<'a> {
    pub interval: Interval,
    pub counts: CountsBufferContent<'a>,
}

#[cfg_attr(test, automock)]
pub trait NucCounter<R: AlignedRead> {
    fn roi(&self) -> &Interval;
    fn process(&mut self, read: &mut R);
    fn reads_counted(&self) -> u32;
    fn reset(&mut self, roi: Interval);
    fn content(&'_ self) -> NucCounterContent<'_>;
    fn empty(&self) -> bool {
        self.reads_counted() == 0
    }
}
