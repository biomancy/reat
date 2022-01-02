use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::mismatches::interval::IntermediateIntervalMismatches;

mod engine;
mod shared;

pub use engine::IntervalStrandingEngine;

pub trait IntervalStrandPredictor<T: IntermediateIntervalMismatches> {
    fn predict(&self, blocks: Vec<T>) -> Vec<T>;
}
