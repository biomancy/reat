use crate::core::hooks::filters::Filter;
use crate::core::hooks::stats::{EditingStat, TypedEditingStat};
use crate::core::hooks::{Hook, HooksEngine};
use crate::core::io::bed::BedRecord;
use crate::core::mismatches::{BatchedMismatches, MismatchesIntermediate};

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

impl<T: BatchedMismatches> Hook<T> for REATHooksEngine<T> {
    fn on_created(&mut self, mismatches: &mut MismatchesIntermediate<T>) {
        for s in &mut self.stats {
            s.on_created(mismatches);
        }
        for s in &mut self.filters {
            s.on_created(mismatches);
        }
    }

    fn on_stranded(&mut self, mismatches: &mut MismatchesIntermediate<T>) {
        for s in &mut self.stats {
            s.on_stranded(mismatches);
        }
        for s in &mut self.filters {
            s.on_stranded(mismatches);
        }
    }

    fn on_finish(&mut self, mismatches: &mut MismatchesIntermediate<T>) {
        for s in &mut self.stats {
            s.on_finish(mismatches);
        }
        for s in &mut self.filters {
            s.on_finish(mismatches);
        }
    }
}

impl<T: BatchedMismatches> HooksEngine<T> for REATHooksEngine<T> {
    fn stats(&self) -> &[Box<dyn EditingStat<T>>] {
        &self.stats
    }

    fn filters(&self) -> &[Box<dyn Filter<T>>] {
        &self.filters
    }
}

// impl<T> HooksEngine<T> for REATHooksEngine<T> {
//     fn add_stat(&mut self, stat: Box<dyn EditingStatHook<T>>) {
//         self.stats.push(stat);
//     }
//
//     fn on_finish(&mut self, mismatches: Vec<T>) -> Vec<T> {
//         // Calculate stats first
//         for h in &mut self.stats {
//             h.hook(&mismatches);
//         }
//
//         let result = if let Some(s) = skip {
//             tofilter.extend(s.into_iter());
//             tofilter
//         } else {
//             tofilter
//         };
//
//         result
//     }
//
//     fn finalize(self) -> Vec<Box<dyn EditingStatHook<T>>> {
//         self.stats
//     }
// }

// pub fn set_skipfilter(&mut self, skipfilter: Box<dyn AlwaysRetainFilter<T>>) {
//     self.skipfilter = Some(skipfilter);
// }
//
// pub fn add_filter(&mut self, filter: Box<dyn Filter<T>>) {
//     self.filters.push(filter);
// }
//
// fn partition(&self, mismatches: Vec<T>) -> (Option<Vec<T>>, Vec<T>) {
//     if let Some(ref s) = self.skipfilter {
//         let (skip, tofilter) = s.partition(mismatches);
//         (Some(skip), tofilter)
//     } else {
//         (None, mismatches)
//     }
// }
