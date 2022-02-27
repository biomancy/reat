#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::hooks::filters::{AlwaysRetainFilter, Filter};
use crate::core::hooks::stats::{EditingStat, EditingStatHook};
use crate::core::io::bed::BedRecord;
use crate::core::mismatches::roi::ROIMismatches;
use crate::core::mismatches::BinnedMismatches;

pub mod engine;
pub mod filters;
pub mod stats;

pub trait Hook<T: BinnedMismatches> {
    fn on_created(&mut self, toretain: &mut Vec<T>, items: &mut Vec<T>);
    fn on_stranded(&mut self, toretain: &mut Vec<T>, items: &mut Vec<T>);
    fn on_finish(&mut self, toretain: &mut Vec<T>, items: &mut Vec<T>);
}

// No master filter for now. All simple filtering will be done prior to storing mismatches (in the counter)
// Retaining now can be defined there as well.
pub trait MasterFilter<T: BinnedMismatches>: Hook<T> {
    type MismatchesPreview;
    fn prefilter(&self, preview: &Self::MismatchesPreview) -> bool;
}

pub trait StatsCalculator<T: BinnedMismatches>: Hook<T> {
    fn results(self) -> ();
}

pub trait HooksEngine<T: BinnedMismatches, F: MasterFilter<T>, S: StatsCalculator<T>>: Hook<T> {
    // Stats are called after filtering at each stage
    fn stats(&self) -> &S;
    fn filter(&self) -> &F;
}
