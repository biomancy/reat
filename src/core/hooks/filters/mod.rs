pub use mismatches::ByMismatches;

use crate::core::mismatches::BatchedMismatches;

use super::Hook;

mod mismatches;

pub trait Filter<T: BatchedMismatches>: Hook<T> {}
