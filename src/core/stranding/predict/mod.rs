use bio_types::strand::Strand;
use dyn_clone::DynClone;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use engine::REATStrandingEngine;

use crate::core::mismatches::{BatchedMismatches, StrandingCounts};

pub use self::ctx::StrandingContext;

pub mod algo;
mod ctx;
mod engine;

pub trait StrandingEngine<T: BatchedMismatches> {
    fn strand(&self, items: Vec<T>) -> Vec<T>;
}

pub enum StrandingAlgoResult {
    AllElements(Strand),
    EachElement((Vec<Strand>, StrandingCounts)),
}

pub trait StrandingAlgo<T: BatchedMismatches>: DynClone + Send {
    fn predict(&self, mismatches: &T) -> StrandingAlgoResult;
}

dyn_clone::clone_trait_object!(<T> StrandingAlgo<T> where T: BatchedMismatches);
