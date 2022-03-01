use itertools::zip;

use crate::core::hooks::filters::Filter;
use crate::core::hooks::Hook;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::roi::{
    BatchedROIMismatches, MismatchesSummary, REATBatchedROIMismatches, REATROIMismatches,
};
use crate::core::mismatches::site::{BinnedSiteMismatches, REATBatchedSiteMismatches};
use crate::core::mismatches::{prefilters, BatchedMismatches, MismatchesIntermediate};

pub struct ByMismatches {
    inner: prefilters::ByMismatches,
}

impl Into<ByMismatches> for prefilters::ByMismatches {
    fn into(self) -> ByMismatches {
        ByMismatches { inner: self }
    }
}

impl ByMismatches {
    fn filter_rois(&self, mm: REATBatchedROIMismatches) -> Option<REATBatchedROIMismatches> {
        let mut oneok = false;
        let mask = mm
            .mismatches()
            .iter()
            .map(|x| {
                let pass = self.inner.ok_roi(x);
                oneok = oneok || pass;
                pass
            })
            .collect();

        match oneok {
            true => Some(mm.filter(mask)),
            false => None,
        }
    }

    fn filter_intervals(&self, mm: REATBatchedSiteMismatches) -> Option<REATBatchedSiteMismatches> {
        let mut oneok = false;
        let mask = zip(mm.prednuc(), mm.seqnuc())
            .map(|(&nuc, &seqnuc)| {
                let pass = self.inner.ok_site(&(nuc, seqnuc));
                oneok = oneok || pass;
                pass
            })
            .collect();

        match oneok {
            true => Some(mm.filter(mask)),
            false => None,
        }
    }
}

impl Hook<REATBatchedROIMismatches> for ByMismatches {
    fn on_stranded(&mut self, mm: &mut MismatchesIntermediate<REATBatchedROIMismatches>) {
        mm.other = std::mem::take(&mut mm.other).into_iter().filter_map(|x| self.filter_rois(x)).collect();
    }
}
impl Filter<REATBatchedROIMismatches> for ByMismatches {}

impl Hook<REATBatchedSiteMismatches> for ByMismatches {
    fn on_stranded(&mut self, mm: &mut MismatchesIntermediate<REATBatchedSiteMismatches>) {
        mm.other = std::mem::take(&mut mm.other).into_iter().filter_map(|x| self.filter_intervals(x)).collect();
    }
}
impl Filter<REATBatchedSiteMismatches> for ByMismatches {}
