mod mismatches;
pub mod retain;

pub use mismatches::ByMismatches;

pub trait MismatchesPreFilter<T> {
    fn is_ok(&self, preview: &T) -> bool;
}
