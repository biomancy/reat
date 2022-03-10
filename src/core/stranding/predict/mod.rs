use bio_types::strand::Strand;
use dyn_clone::DynClone;

pub use engine::REATStrandingEngine;

use crate::core::mismatches::{MismatchesVec, StrandingCounts};

use crate::core::strandutil::StrandedData;

pub mod algo;
mod engine;

pub trait StrandingEngine<T: MismatchesVec> {
    fn strand(&self, items: StrandedData<Option<T>>) -> StrandedData<Option<T>>;
}

pub enum StrandingAlgoResult {
    AllElements(Strand),
    EachElement((Vec<Strand>, StrandingCounts)),
}

pub trait StrandingAlgo<T: MismatchesVec>: DynClone + Send {
    fn predict(&self, mismatches: &T) -> StrandingAlgoResult;
}

dyn_clone::clone_trait_object!(<T> StrandingAlgo<T> where T: MismatchesVec);
