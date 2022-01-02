mod interval;
mod roi;

pub use interval::{IntervalStrandPredictor, IntervalStrandingEngine};
pub use roi::{ROIStrandPredictor, ROIStrandingEngine};

use bio_types::strand::Strand;
#[cfg(test)]
use mockall::{automock, predicate::*};
