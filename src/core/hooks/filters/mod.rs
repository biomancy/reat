pub use mismatches::ByMismatches;

use crate::core::mismatches::MismatchesVec;

use super::Hook;

mod mismatches;
mod utils;

pub trait Filter<T: MismatchesVec>: Hook<T> {}
