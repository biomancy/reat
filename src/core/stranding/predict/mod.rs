#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::mismatches::roi::IntermediateROIMismatches;
use crate::core::mismatches::IntermediateMismatches;

pub mod interval;
pub mod roi;
pub mod shared;

// TODO: refactor once we have stable trait specialization
pub trait StrandingEngine<T: IntermediateMismatches> {
    // TODO: It's a lot better to pass CONST objects but SiteBlocks can be splitted during stranding
    fn strand(&self, items: Vec<T>) -> Vec<T>;
}
