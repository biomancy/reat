use bio_types::genome::Interval;

pub use collider::ReadsCollider;

use crate::core::read::AlignedRead;

mod collider;
pub mod hts;
pub mod ncounters;

// #[cfg_attr(test, automock)]
// Pileup engine
pub trait ReadsPileupEngine<'a, R: AlignedRead, Collider: ReadsCollider<'a, R>> {
    // Reset and run the engine and collider for the given interval and get results
    fn run(&mut self, interval: Interval, cwork: Collider::Workload);
    // TODO: Proper error handling
    // A separate method handler with different const self

    // TODO: return not a vec, but simply a collider result;
    // Didn't find a way to override this type
    // TODO: create an issue
    fn result(&'a self) -> Result<Collider::ColliderResult, ()>;
}
