use bio_types::genome::Interval;

pub use base::BaseRunCtx;

use crate::rada::counting::NucCounterContent;
use crate::rada::dna::Nucleotide;
use crate::rada::filtering::summary::{IntervalSummaryFilter, LocusSummaryFilter};
use crate::rada::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};
use std::path::PathBuf;

mod base;

pub trait RunCtx {
    fn count(&mut self, interval: &Interval);
    fn finalize(&self) -> Option<(NucCounterContent, Vec<Nucleotide>)>;
    fn htsfiles(&self) -> &[PathBuf];
}

pub trait ROIRunCtx: RunCtx {
    type StrandPred: IntervalStrandPredictor;
    type Filter: IntervalSummaryFilter;

    fn strandpred(&self) -> &Self::StrandPred;
    fn filter(&self) -> &Self::Filter;
}

pub trait LociRunCtx: RunCtx {
    type StrandPred: LocusStrandPredictor;
    type Filter: LocusSummaryFilter;

    fn strandpred(&self) -> &Self::StrandPred;
    fn filter(&self) -> &Self::Filter;
}
