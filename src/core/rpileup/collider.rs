use crate::core::read::AlignedRead;
use bio_types::genome::AbstractInterval;

// #[cfg_attr(test, automock)]
// A function computed on top of sequenced filters in a given interval
pub trait ReadsCollider<'a, R: AlignedRead> {
    type ColliderResult;
    type Workload: AbstractInterval;

    fn reset(&mut self, info: Self::Workload);
    fn collide(&mut self, read: &R);
    fn finalize(&mut self);
    fn result(&'a self) -> Self::ColliderResult;
}
