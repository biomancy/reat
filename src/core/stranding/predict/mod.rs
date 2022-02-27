use derive_getters::Dissolve;
#[cfg(test)]
use mockall::{automock, predicate::*};

pub use engine::REATStrandingEngine;

use crate::core::mismatches::roi::BinnedROIMismatches;
use crate::core::mismatches::site::BinnedSiteMismatches;
use crate::core::mismatches::BatchedMismatches;

pub use self::ctx::StrandingContext;

pub mod algo;
mod ctx;
mod engine;

pub trait StrandingEngine<T: BatchedMismatches> {
    fn add_algo(&mut self, algo: Box<dyn StrandingAlgo<T>>);
    fn strand(&self, items: Vec<T>) -> Vec<T>;
}

pub trait StrandingAlgo<T: BatchedMismatches> {
    fn predict(&self, ctx: StrandingContext<T>) -> StrandingContext<T>;
}
