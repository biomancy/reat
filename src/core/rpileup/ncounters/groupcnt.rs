use std::marker::PhantomData;

use bio_types::genome::{AbstractInterval, Interval};

use super::CountingResults;
use derive_getters::Getters;

#[derive(Getters)]
pub struct GroupedNucCounts<'a, T: CountingResults<'a>> {
    interval: Interval,
    items: Vec<T>,
    marker: PhantomData<&'a ()>,
}

impl<'a, T: CountingResults<'a>> GroupedNucCounts<'a, T> {
    pub fn new(interval: Interval, items: Vec<T>) -> Self {
        debug_assert!(items.iter().all(|x| x.contig() == interval.contig()
            && interval.range().contains(&x.range().start)
            && interval.range().contains(&x.range().end)));
        Self { interval, items, marker: Default::default() }
    }

    pub fn dissolve(self) -> (Interval, Vec<T>) {
        (self.interval, self.items)
    }
}
