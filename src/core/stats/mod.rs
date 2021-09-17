#[cfg(test)]
use mockall::mock;

pub use editing_index::EditingIndex;

use crate::core::summary::ROISummary;

mod editing_index;

pub trait ROIBasedStat: Sized {
    fn process(&mut self, summary: &ROISummary);
    fn combine(stats: &[Self]) -> Self;
    // TODO in the future it could be a trait (e.g. TableDisplay)
    fn header() -> &'static str;
    fn row(&self) -> String;
}

#[cfg(test)]
mock! {
    pub ROIBasedStat {}
    impl ROIBasedStat for ROIBasedStat {
        fn process(&mut self, summary: &ROISummary);
        fn combine(stats: &[MockROIBasedStat]) -> MockROIBasedStat;
        fn header() -> &'static str;
        fn row(&self) -> String;
    }
}
