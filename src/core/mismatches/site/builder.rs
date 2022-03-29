use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;
use itertools::izip;

use crate::core::dna::NucCounts;
use crate::core::mismatches::prefilters::retain::SitesRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::site::{SiteData, SiteDataVec, SiteMismatchesVec};
use crate::core::mismatches::Batch;
use crate::core::refpred::{RefEngine, RefEngineResult};
use crate::core::rpileup::ncounter::NucCounterResult;
use crate::core::strandutil::Stranded;

use super::super::Builder;

#[derive(Clone)]
pub struct SiteMismatchesBuilder<SR, MP> {
    buffer: Vec<NucCounts>,
    refpred: Box<dyn RefEngine>,
    retainer: Option<SR>,
    prefilter: Option<MP>,
}

impl<'a, SR, MP> SiteMismatchesBuilder<SR, MP>
where
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteData>,
{
    pub fn new(maxsize: usize, refpred: Box<dyn RefEngine>, retainer: Option<SR>, prefilter: Option<MP>) -> Self {
        Self { buffer: Vec::with_capacity(maxsize), refpred, retainer, prefilter }
    }

    fn process(
        &self,
        retained: &[Range<Position>],
        cntrange: Range<Position>,
        cnts: &[NucCounts],
        refngn: &RefEngineResult,
        retbuilder: &mut SiteDataVec,
        othbuilder: &mut SiteDataVec,
    ) {
        debug_assert_eq!(cntrange.end - cntrange.start, cnts.len() as Position);
        debug_assert_eq!(cnts.len(), refngn.reference.len());
        debug_assert_eq!(cnts.len(), refngn.predicted.len());
        debug_assert!(retained.iter().all(|x| cntrange.contains(&x.start) && cntrange.contains(&x.end)));

        let retsize = retained.iter().map(|x| x.end - x.start).sum::<Position>() as usize;
        debug_assert!(retsize <= cnts.len());

        // Retain iterator
        let mut reiter = retained.iter();
        let mut retrange = reiter.next();

        for (pos, (&cnt, &refnuc, &prednuc)) in izip!(cnts, refngn.reference, refngn.predicted).enumerate() {
            let pos = pos as Position + cntrange.start;
            // Do we need to move the iterator?
            if retrange.map_or(false, |x| x.end <= pos) {
                retrange = reiter.next();
            }

            // Are we inside the retained region?
            let data = SiteData { pos, refnuc, prednuc, sequenced: cnt };
            if retrange.map_or(false, |x| x.contains(&pos)) {
                retbuilder.push(data);
            } else if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&data)) {
                othbuilder.push(data);
            }
        }
    }

    #[inline]
    fn size_hint(&self, nc: &NucCounterResult<'a, ()>) -> Stranded<usize> {
        let mut size: Stranded<usize> = Stranded::default();
        for item in nc.cnts.iter() {
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                size[strand] += item.cnts[strand].map_or(0, |x| x.len());
            }
        }
        size
    }
}

impl<'a, SR, MP> Builder<'a> for SiteMismatchesBuilder<SR, MP>
where
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteData>,
{
    type Out = SiteMismatchesVec;
    type SourceCounts = NucCounterResult<'a, ()>;

    fn build(&mut self, nc: Self::SourceCounts) -> Batch<Self::Out> {
        let contig = nc.contig;

        // Pre-allocate results
        let _hint = self.size_hint(&nc);
        let mut items = Stranded::with_fn(|strnd| {
            // SiteMismatchesVec::new(contig.to_owned(), strnd, SiteDataVec::with_capacity(hint[strnd]))
            SiteMismatchesVec::new(contig.to_owned(), strnd, SiteDataVec::new())
        });
        let mut retained = Stranded::with_fn(|strnd| {
            SiteMismatchesVec::new(contig.to_owned(), strnd, SiteDataVec::new())
            // SiteMismatchesVec::new(contig.to_owned(), strnd, SiteDataVec::with_capacity(hint[strnd] / 10))
        });

        for item in nc.cnts.into_iter() {
            // Predict the reference
            let counts = item.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(contig, item.range.clone(), counts);
            let reference = self.refpred.results();

            // Find loci that must be retained
            let mustloci = self.retainer.as_ref().map_or(vec![], |r| r.retained(contig, item.range.clone()));
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                if let Some(cnt) = item.cnts[strand] {
                    // debug_assert!(item.coverage[strand] > 0);
                    self.process(
                        &mustloci,
                        item.range.clone(),
                        cnt,
                        &reference,
                        &mut retained[strand].data,
                        &mut items[strand].data,
                    );
                };
            }
        }

        Batch { contig: contig.to_owned(), mapped: nc.mapped, retained, items }
    }
}
