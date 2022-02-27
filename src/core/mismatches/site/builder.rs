use crate::core::mismatches::prefilters::retain::SitesRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::site::REATBatchedSiteMismatches;
use crate::core::mismatches::PreFilteredBatchedMismatches;
use crate::core::refpred::RefEngine;
use crate::core::rpileup::ncounters::cnt::IntervalNucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCounts;

use super::super::BatchedMismatchesBuilder;
use super::SiteMismatchesPreview;
use crate::core::dna::NucCounts;
use bio_types::genome::AbstractInterval;
use itertools::zip;

pub struct REATSiteMismatchesBuilder<RE: RefEngine, SR: SitesRetainer, MP: MismatchesPreFilter<SiteMismatchesPreview>> {
    buffer: Vec<NucCounts>,
    refpred: RE,
    retainer: Option<SR>,
    prefilter: Option<MP>,
}

impl<'a, RE, SR, MP> BatchedMismatchesBuilder<'a, IntervalNucCounts<'a>> for REATSiteMismatchesBuilder<RE, SR, MP>
where
    RE: RefEngine,
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteMismatchesPreview>,
{
    type Mismatches = REATBatchedSiteMismatches;

    fn build(&mut self, nc: IntervalNucCounts<'a>) -> PreFilteredBatchedMismatches<Self::Mismatches> {
        let (retained, other) = (Vec::new(), Vec::new());
        let contig = nc.contig();

        for nc in nc.consume() {
            debug_assert!(nc.forward.is_some() || nc.reverse.is_some() || nc.unstranded.is_some());
            // Gather sequenced nucleotides in each position
            let counts = match (nc.forward, nc.reverse, nc.unstranded) {
                (Some(c), None, None) => c,
                (None, Some(c), None) => c,
                (None, None, Some(c)) => c,
                (_, _, _) => {
                    self.buffer.clear();
                    self.buffer.resize(nc.range.end as usize - nc.range.start as usize, Default::default());
                    for x in [nc.forward, nc.reverse, nc.unstranded] {
                        if let Some(x) = x {
                            debug_assert!(x.len() == self.buffer.len());
                            for (bc, xc) in zip(self.buffer.iter_mut(), x.iter()) {
                                *bc += *xc;
                            }
                        }
                    }
                    &self.buffer
                }
            };

            // Predict the reference
            self.refpred.reset();
            self.refpred.run(contig, nc.range, counts);
            let reference = self.refpred.results();

            // Special case -> unstranded reads
        }

        PreFilteredBatchedMismatches { retained, other }
    }
}
