use std::ops::Add;

#[cfg(test)]
use mockall::mock;

pub use editing_index::EditingIndex;

use crate::rada::summary::IntervalSummary;

mod editing_index;

pub trait IntervalBasedStat: Add<Output = Self> + Default {
    fn process(&mut self, summary: &IntervalSummary);

    // TODO in the future it could be something like a trait
    fn header() -> &'static str;
    fn row(&self) -> String;
}

#[cfg(test)]
mock! {
    pub IntervalBasedStat {}

    impl Add<MockIntervalBasedStat> for IntervalBasedStat {
        type Output = MockIntervalBasedStat;

        fn add(self, rhs: MockIntervalBasedStat) -> MockIntervalBasedStat;
    }

    impl IntervalBasedStat for IntervalBasedStat {
        fn process(&mut self, summary: &IntervalSummary);

        fn header() -> &'static str;
        fn row(&self) -> String;
    }
}
