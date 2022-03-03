use std::any::Any;

#[cfg(test)]
use mockall::mock;

pub use roi_editing_index::ROIEditingIndex;

use crate::core::hooks::Hook;
use crate::core::mismatches::BatchedMismatches;

mod roi_editing_index;

#[derive(Hash, PartialEq, Eq)]
pub enum EditingStatType {
    ROIEditingIndex,
}

pub trait EditingStat<T: BatchedMismatches>: Hook<T> + Any {
    fn into_any(self: Box<Self>) -> (EditingStatType, Box<dyn Any>);
}
