use bio_types::strand::Strand;
use dyn_clone::DynClone;

pub use engine::REATStrandingEngine;

use crate::core::mismatches::{MismatchesVec, StrandingCounts};

use crate::core::strandutil::Stranded;

pub mod algo;
mod engine;

pub trait StrandingEngine<T: MismatchesVec> {
    fn strand(&self, contig: &str, items: Stranded<T>) -> Stranded<T>;
}

pub trait StrandingAlgo<T: MismatchesVec>: DynClone + Send {
    fn predict(&self, contig: &str, items: &mut Stranded<T>);
}

dyn_clone::clone_trait_object!(<T> StrandingAlgo<T> where T: MismatchesVec);
