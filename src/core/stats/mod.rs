use std::ops::Add;

#[cfg(test)]
use mockall::automock;

pub use editing_index::EditingIndex;

use crate::core::summary::ROISummary;

mod editing_index;

// #[cfg_attr(test, automock)]
pub trait ROIBasedStat: Sized {
    fn process(&mut self, summary: &ROISummary);
    fn combine(stats: &[Self]) -> Self;
    // TODO in the future it could be a trait (e.g. TableDisplay)
    fn header() -> &'static str;
    fn row(&self) -> String;
}
