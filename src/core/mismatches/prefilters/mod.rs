pub use mismatches::ByMismatches;

mod mismatches;
pub mod retain;

pub trait MismatchesPreFilter<T> {
    fn is_ok(&self, preview: &T) -> bool;
}
