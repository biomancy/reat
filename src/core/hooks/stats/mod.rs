use std::fmt::Display;

use crate::core::hooks::Hook;
use crate::core::mismatches::BatchedMismatches;
#[cfg(test)]
use mockall::mock;

pub use roi_editing_index::ROIEditingIndex;

mod roi_editing_index;

pub enum TypedEditingStat {
    ROIEditingIndex(Box<ROIEditingIndex>),
}

pub trait EditingStat<T: BatchedMismatches>: Hook<T> {
    fn downcast(self: Box<Self>) -> TypedEditingStat;
}
