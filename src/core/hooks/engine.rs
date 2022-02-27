use crate::core::hooks::filters::{AlwaysRetainFilter, Filter};
use crate::core::hooks::stats::{EditingStat, EditingStatHook};
use crate::core::hooks::HooksEngine;
use crate::core::io::bed::BedRecord;

pub struct REATHooksEngine<T> {
    stats: Vec<Box<dyn EditingStatHook<T>>>,
}

impl<T> REATHooksEngine<T> {
    pub fn new() -> Self {
        Self { stats: vec![], skipfilter: None, filters: vec![] }
    }
}

impl<T> HooksEngine<T> for REATHooksEngine<T> {
    fn add_stat(&mut self, stat: Box<dyn EditingStatHook<T>>) {
        self.stats.push(stat);
    }

    fn on_finish(&mut self, mismatches: Vec<T>) -> Vec<T> {
        // Calculate stats first
        for h in &mut self.stats {
            h.hook(&mismatches);
        }

        let result = if let Some(s) = skip {
            tofilter.extend(s.into_iter());
            tofilter
        } else {
            tofilter
        };

        result
    }

    fn finalize(self) -> Vec<Box<dyn EditingStatHook<T>>> {
        self.stats
    }
}
