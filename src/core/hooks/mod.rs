use dyn_clone::DynClone;
#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::hooks::stats::EditingStat;
use crate::core::mismatches::{BatchedMismatches, FilteredBatchedMismatches};

pub mod engine;
pub mod filters;
pub mod stats;

pub trait Hook<T: BatchedMismatches>: DynClone + Send {
    fn on_finish(&mut self, _mismatches: &mut FilteredBatchedMismatches<T>) {}
}

dyn_clone::clone_trait_object!(<T> Hook<T> where T: BatchedMismatches);

// Stats must be called after filtering at each stage
pub trait HooksEngine<T: BatchedMismatches>: Hook<T> {
    fn stats(self) -> Vec<Box<dyn EditingStat<T>>>;
}
