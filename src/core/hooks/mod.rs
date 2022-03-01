#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::hooks::filters::Filter;
use crate::core::hooks::stats::{EditingStat, TypedEditingStat};
use crate::core::mismatches::roi::ROIMismatches;
use crate::core::mismatches::{BatchedMismatches, MismatchesIntermediate};

pub mod engine;
pub mod filters;
pub mod stats;

pub trait Hook<T: BatchedMismatches> {
    fn on_created(&mut self, mismatches: &mut MismatchesIntermediate<T>) {}
    fn on_stranded(&mut self, mismatches: &mut MismatchesIntermediate<T>) {}
    fn on_finish(&mut self, mismatches: &mut MismatchesIntermediate<T>) {}
}

// Stats must be called after filtering at each stage
pub trait HooksEngine<T: BatchedMismatches>: Hook<T> {
    fn stats(&self) -> &[Box<dyn EditingStat<T>>];
    fn filters(&self) -> &[Box<dyn Filter<T>>];
}
