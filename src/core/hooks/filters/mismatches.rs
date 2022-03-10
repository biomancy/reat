use itertools::zip;

use crate::core::hooks::filters::{utils, Filter};
use crate::core::hooks::Hook;
use crate::core::mismatches::roi::{REATROIMismatchesVec, ROIMismatches, ROIMismatchesVec};
use crate::core::mismatches::site::{REATSiteMismatchesVec, SiteMismatches, SiteMismatchesVec};
use crate::core::mismatches::{prefilters, Context, MismatchesVec};

#[derive(Clone)]
pub struct ByMismatches {
    inner: prefilters::ByMismatches,
}

impl Into<ByMismatches> for prefilters::ByMismatches {
    fn into(self) -> ByMismatches {
        ByMismatches { inner: self }
    }
}

impl ByMismatches {
    fn filter_rois<T>(&self, mm: T) -> Option<T>
    where
        T: ROIMismatchesVec,
        T::Flat: ROIMismatches,
    {
        // let iter = mm.mismatches().iter().map(|x| self.inner.enough_mismatches_per_roi(x));
        // utils::maskiter(iter).map(|x| mm.filter(x))
        todo!()
    }

    fn filter_intervals<T>(&self, mm: T) -> Option<T>
    where
        T: SiteMismatchesVec,
        T::Flat: SiteMismatches,
    {
        // let iter =
        //     zip(mm.prednuc(), mm.seqnuc()).map(|(&nuc, &seqnuc)| self.inner.enough_mismatches_per_site(&(nuc, seqnuc)));
        // utils::maskiter(iter).map(|x| mm.filter(x))
        todo!()
    }
}

impl Hook<REATROIMismatchesVec> for ByMismatches {
    fn on_finish(&mut self, mm: &mut Context<REATROIMismatchesVec>) {
        mm.items.apply_mut(|opt, _| opt.and_then(|x| self.filter_rois(x)));
    }
}

impl Filter<REATROIMismatchesVec> for ByMismatches {}

impl Hook<REATSiteMismatchesVec> for ByMismatches {
    fn on_finish(&mut self, mm: &mut Context<REATSiteMismatchesVec>) {
        mm.items.apply_mut(|opt, _| opt.and_then(|x| self.filter_intervals(x)));
    }
}

impl Filter<REATSiteMismatchesVec> for ByMismatches {}
