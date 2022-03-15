use bio_types::genome::AbstractInterval;

use crate::core::read::AlignedRead;

pub mod hts;
pub mod ncounter;

// Pileup engine
pub trait ReadsCollidingEngine<R: AlignedRead, Collider>
where
    R: AlignedRead,
    Collider: for<'a> ReadsCollider<'a, R>,
{
    // Reset and run the engine and collider for the given interval and get results
    fn run(&mut self, cwork: <Collider as ReadsCollider<'_, R>>::Workload);
    // Get calculated result if any is available
    fn result(&self) -> Option<<Collider as ReadsCollider<'_, R>>::ColliderResult>;
}

// A function computed on top of sequenced filters in a given interval
pub trait ReadsCollider<'a, R: AlignedRead> {
    type ColliderResult;
    type Workload: AbstractInterval;

    // Reset the collider using the given Workload
    fn reset(&mut self, info: Self::Workload);
    // Run the collider
    fn collide(&mut self, read: &R);
    // Calculate the result
    fn finalize(&mut self);
    // Return prepared info
    fn result(&'a self) -> Self::ColliderResult;
}
