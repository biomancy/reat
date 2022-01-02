use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{automock, predicate::*};

use crate::core::mismatches::roi::IntermediateROIMismatches;

mod engine;
mod shared;
pub use engine::ROIStrandingEngine;

pub trait ROIStrandPredictor<T: IntermediateROIMismatches> {
    fn predict(&self, rois: Vec<T>) -> Vec<T>;
}
