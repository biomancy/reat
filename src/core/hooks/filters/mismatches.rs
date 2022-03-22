use crate::core::hooks::filters::Filter;
use crate::core::hooks::Hook;
use crate::core::mismatches::roi::ROIMismatchesVec;
use crate::core::mismatches::site::SiteMismatchesVec;
use crate::core::mismatches::{prefilters, Batch};
use crate::core::refpred::PredNucleotide;

#[derive(Clone)]
pub struct ByMismatches {
    inner: prefilters::ByMismatches,
}

impl From<prefilters::ByMismatches> for ByMismatches {
    fn from(bm: prefilters::ByMismatches) -> Self {
        Self { inner: bm }
    }
}

impl Hook<ROIMismatchesVec> for ByMismatches {
    fn on_finish(&mut self, mm: &mut Batch<ROIMismatchesVec>) {
        mm.items.apply_mut(|x, _| x.data.retain(|x| self.inner.enough_mismatches_per_roi(x.mismatches)));
    }
}

impl Filter<ROIMismatchesVec> for ByMismatches {}

impl Hook<SiteMismatchesVec> for ByMismatches {
    fn on_finish(&mut self, mm: &mut Batch<SiteMismatchesVec>) {
        mm.items.apply_mut(|x, _| {
            x.data.retain(|x| match x.prednuc {
                PredNucleotide::Homozygous(nuc) => self.inner.enough_mismatches_per_site(*nuc, x.sequenced),
                PredNucleotide::Heterozygous((n1, n2)) => {
                    self.inner.enough_mismatches_per_site(*n1, x.sequenced)
                        || self.inner.enough_mismatches_per_site(*n2, x.sequenced)
                }
            })
        });
    }
}

impl Filter<SiteMismatchesVec> for ByMismatches {}
