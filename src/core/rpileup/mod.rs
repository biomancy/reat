use bio_types::genome::Interval;

pub use collider::ReadsCollider;

use crate::core::read::AlignedRead;

mod collider;
pub mod hts;
pub mod ncounters;

// #[cfg_attr(test, automock)]
// Pileup engine
pub trait ReadsCollidingEngine<'a, R: AlignedRead, Collider: ReadsCollider<'a, R>> {
    // Reset and run the engine and collider for the given interval and get results
    fn run(&mut self, cwork: Collider::Workload);
    // TODO: Proper error handling
    fn result(&'a self) -> Result<Collider::ColliderResult, ()>;
}
