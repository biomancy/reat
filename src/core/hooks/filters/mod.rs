pub use mismatches::ByMismatches;

use super::Hook;
use crate::core::mismatches::BatchedMismatches;

mod mismatches;

pub trait Filter<T: BatchedMismatches>: Hook<T> {}
