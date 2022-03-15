use std::any::Any;

pub use roi_editing_index::ROIEditingIndex;

use crate::core::hooks::Hook;
use crate::core::mismatches::MismatchesVec;

mod roi_editing_index;

#[derive(Hash, PartialEq, Eq)]
pub enum EditingStatType {
    ROIEditingIndex,
}

pub trait EditingStat<T: MismatchesVec>: Hook<T> + Any {
    fn into_any(self: Box<Self>) -> (EditingStatType, Box<dyn Any>);
}
