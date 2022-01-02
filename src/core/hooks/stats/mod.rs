use std::fmt::Display;

#[cfg(test)]
use mockall::mock;

pub use roi_editing_index::ROIEditingIndex;

pub use crate::core::hooks::ProcessingHook;

mod roi_editing_index;

pub trait EditingStat<T>: ProcessingHook<T> {
    fn combine(stats: &[Self]) -> Self
    where
        Self: Sized;
}

// #[cfg(test)]
// mock! {
//     pub ROIBasedStat {}
//     impl ROIBasedStat for ROIBasedStat {
//         fn process(&mut self, mismatches: &ROISummary);
//         fn combine(stats: &[MockROIBasedStat]) -> MockROIBasedStat;
//         fn header() -> &'static str;
//         fn row(&self) -> String;
//     }
// }
