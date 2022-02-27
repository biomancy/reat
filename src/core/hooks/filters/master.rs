use crate::core::hooks::{Filter, Hook, MasterFilter};
use crate::core::mismatches::BinnedMismatches;

use super::AlwaysRetainFilter;

pub struct REATMasterFilter<T: BinnedMismatches> {
    skipfilter: Option<Box<dyn AlwaysRetainFilter<T>>>,
    filters: Vec<Box<dyn Filter<T>>>,
}

impl<T: BinnedMismatches> REATMasterFilter<T> {
    pub fn set_skipfilter(&mut self, skipfilter: Box<dyn AlwaysRetainFilter<T>>) {
        self.skipfilter = Some(skipfilter);
    }

    pub fn add_filter(&mut self, filter: Box<dyn Filter<T>>) {
        self.filters.push(filter);
    }

    fn partition(&self, mismatches: Vec<T>) -> (Option<Vec<T>>, Vec<T>) {
        if let Some(ref s) = self.skipfilter {
            let (skip, tofilter) = s.partition(mismatches);
            (Some(skip), tofilter)
        } else {
            (None, mismatches)
        }
    }
}

impl<T: BinnedMismatches> Hook<T> for REATMasterFilter<T> {
    fn on_created(&mut self, mismatches: Vec<T>) -> Vec<T> {
        if self.filters.is_empty() {
            return mismatches;
        }
        // Separate force elements
        let (skip, mut tofilter) = self.partition(mismatches);

        for h in &mut self.filters {
            tofilter = h.hook(tofilter);
        }
    }

    fn on_stranded(&mut self, mismatches: Vec<T>) -> Vec<T> {
        if self.filters.is_empty() {
            return mismatches;
        }
        todo!()
    }

    fn on_finish(&mut self, mismatches: Vec<T>) -> Vec<T> {
        if self.filters.is_empty() {
            return mismatches;
        }
        todo!()
    }
}

impl<T: BinnedMismatches> MasterFilter<T> for REATMasterFilter<T> {
    type MismatchesPreview = ();

    fn prefilter(&self, preview: &Self::MismatchesPreview) -> bool {
        todo!()
    }
}
