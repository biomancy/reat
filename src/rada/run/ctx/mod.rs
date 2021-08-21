use std::path::PathBuf;

use bio_types::genome::Interval;
#[cfg(test)]
use mockall::mock;

pub use base::BaseRunCtx;

use crate::rada::counting::NucCounterContent;
use crate::rada::dna::Nucleotide;
use crate::rada::filtering::summary::{IntervalSummaryFilter, LocusSummaryFilter};
#[cfg(test)]
use crate::rada::filtering::summary::{MockIntervalSummaryFilter, MockLocusSummaryFilter};
use crate::rada::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};
#[cfg(test)]
use crate::rada::stranding::predict::{MockIntervalStrandPredictor, MockLocusStrandPredictor};

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

#[cfg(test)]
mock! {
    pub ROIRunCtx {}

    impl RunCtx for ROIRunCtx {
        fn count(&mut self, interval: &Interval);
        fn finalize<'a>(&'a self) -> Option<(NucCounterContent<'a>, Vec<Nucleotide>)>;
        fn htsfiles(&self) -> &[PathBuf];
    }

    impl ROIRunCtx for ROIRunCtx {
        type StrandPred = MockIntervalStrandPredictor;
        type Filter = MockIntervalSummaryFilter;

        fn strandpred(&self) -> &MockIntervalStrandPredictor;
        fn filter(&self) -> &MockIntervalSummaryFilter;
    }
}

#[cfg(test)]
mock! {
    pub LociRunCtx {}

    impl RunCtx for LociRunCtx {
        fn count(&mut self, interval: &Interval);
        fn finalize<'a>(&'a self) -> Option<(NucCounterContent<'a>, Vec<Nucleotide>)>;
        fn htsfiles(&self) -> &[PathBuf];
    }

    impl LociRunCtx for LociRunCtx {
        type StrandPred = MockLocusStrandPredictor;
        type Filter = MockLocusSummaryFilter;

        fn strandpred(&self) -> &MockLocusStrandPredictor;
        fn filter(&self) -> &MockLocusSummaryFilter;
    }
}
