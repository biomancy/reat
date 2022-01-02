#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::mismatches::interval::IntervalMismatches;
use crate::core::mismatches::roi::ROIMismatches;

pub mod filters;
pub mod stats;

#[cfg_attr(test, automock(type T=();))]
pub trait ProcessingHook<T> {
    fn hook(&mut self, objects: Vec<T>) -> Vec<T>;
}
