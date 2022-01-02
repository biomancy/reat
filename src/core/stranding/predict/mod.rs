#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::mismatches::roi::IntermediateROIMismatches;
use crate::core::mismatches::IntermediateMismatches;

pub mod algo;
pub mod engines;

pub trait StrandingEngine<T: IntermediateMismatches> {
    // TODO: It's a lot better to pass const objects but blocked mismatches can be split during stranding
    fn strand(&self, items: Vec<T>) -> Vec<T>;
}
