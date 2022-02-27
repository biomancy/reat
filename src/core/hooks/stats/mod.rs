use std::fmt::Display;

#[cfg(test)]
use mockall::mock;

pub use roi_editing_index::ROIEditingIndex;

mod roi_editing_index;

pub enum EditingStat {
    ROIEditingIndex,
}

pub trait EditingStatHook<T> {
    fn stype(&self) -> EditingStat;
    fn hook(&mut self, objects: &[T]);
}
