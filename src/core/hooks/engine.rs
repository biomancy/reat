use crate::core::hooks::{Hook, HooksEngine};
use crate::core::hooks::filters::Filter;
use crate::core::hooks::stats::EditingStat;
use crate::core::mismatches::{BatchedMismatches, FilteredBatchedMismatches};

pub struct REATHooksEngine<T> {
    stats: Vec<Box<dyn EditingStat<T>>>,
    filters: Vec<Box<dyn Filter<T>>>,
}

impl<T> REATHooksEngine<T> {
    pub fn new() -> Self {
        Self { stats: vec![], filters: vec![] }
    }

    pub fn add_stat(&mut self, stat: Box<dyn EditingStat<T>>) {
        self.stats.push(stat);
    }

    pub fn add_filter(&mut self, filter: Box<dyn Filter<T>>) {
        self.filters.push(filter);
    }
}

impl<T> Clone for REATHooksEngine<T> {
    fn clone(&self) -> Self {
        Self {
            stats: self.stats.iter().map(|x| dyn_clone::clone_box(x.as_ref())).collect(),
            filters: self.filters.iter().map(|x| dyn_clone::clone_box(x.as_ref())).collect(),
        }
    }
}

impl<T: BatchedMismatches> Hook<T> for REATHooksEngine<T> {
    fn on_finish(&mut self, mismatches: &mut FilteredBatchedMismatches<T>) {
        for s in &mut self.stats {
            s.on_finish(mismatches);
        }
        for f in &mut self.filters {
            f.on_finish(mismatches);
        }
    }
}

impl<T: BatchedMismatches> HooksEngine<T> for REATHooksEngine<T> {
    fn stats(self) -> Vec<Box<dyn EditingStat<T>>> {
        self.stats
    }
}
