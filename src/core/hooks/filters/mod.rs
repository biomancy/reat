pub use mismatches::ByMismatches;

use super::Hook;

mod master;
mod mismatches;
mod retain;

pub trait Filter<T>: Hook<T> {}

pub trait AlwaysRetainFilter<T> {
    fn partition(&self, objects: Vec<T>) -> (Vec<T>, Vec<T>);
}
