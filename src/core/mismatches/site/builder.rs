use std::ops::Range;

use bio_types::genome::AbstractInterval;
use bio_types::genome::Position;
use bio_types::strand::Strand;
use itertools::izip;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::SitesRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::site::REATBatchedSiteMismatches;
use crate::core::mismatches::MismatchesIntermediate;
use crate::core::refpred::RefEngine;
use crate::core::rpileup::ncounters::cnt::IntervalNucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCounts;

use super::super::BatchedMismatchesBuilder;
use super::SiteMismatchesPreview;

pub struct REATSiteMismatchesBuilder<RE: RefEngine, SR: SitesRetainer, MP: MismatchesPreFilter<SiteMismatchesPreview>> {
    buffer: Vec<NucCounts>,
    refpred: RE,
    retainer: Option<SR>,
    prefilter: Option<MP>,
}

impl<'a, RE, SR, MP> REATSiteMismatchesBuilder<RE, SR, MP>
where
    RE: RefEngine,
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteMismatchesPreview>,
{
    fn process(
        &self,
        contig: &str,
        range: Range<Position>,
        strand: Strand,
        retained: &[Range<Position>],
        counts: &[NucCounts],
        reference: &[Nucleotide],
        predref: &[Nucleotide],
    ) -> (Option<REATBatchedSiteMismatches>, Option<REATBatchedSiteMismatches>) {
        let retsize = retained.iter().map(|x| x.end - x.start).sum::<Position>() as usize;
        debug_assert!(retsize <= counts.len());
        let (mut retbuilder, mut othbuilder) =
            (HelperBuilder::new(retsize), HelperBuilder::new(counts.len() - retsize));

        // Retain iterator
        let mut reiter = retained.iter();
        let mut retrange = reiter.next();

        for (loc, (&cnt, &refn, &predn)) in izip!(counts, reference, predref).enumerate() {
            let loc = loc as Position + range.start;
            // Do we need to move the iterator?
            if !retrange.map_or(true, |x| x.contains(&loc)) {
                retrange = reiter.next();
            }

            // Are we inside the retained region?
            if retrange.map_or(false, |x| x.contains(&loc)) {
                retbuilder.add(loc, refn, predn, cnt);
            } else {
                let preview = (predn, cnt);
                if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&preview)) {
                    othbuilder.add(loc, refn, preview.0, preview.1);
                }
            }
        }

        (retbuilder.finalize(contig.into(), strand), othbuilder.finalize(contig.into(), strand))
    }
}

impl<'a, RE, SR, MP> BatchedMismatchesBuilder<'a, IntervalNucCounts<'a>> for REATSiteMismatchesBuilder<RE, SR, MP>
where
    RE: RefEngine,
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteMismatchesPreview>,
{
    type Mismatches = REATBatchedSiteMismatches;

    fn build(&mut self, nc: IntervalNucCounts<'a>) -> MismatchesIntermediate<Self::Mismatches> {
        let (mut retained, mut other) = (Vec::new(), Vec::new());

        let (contig, items) = nc.consume();
        for nc in items {
            debug_assert!(nc.forward.is_some() || nc.reverse.is_some() || nc.unstranded.is_some());

            // Predict the reference
            let counts = nc.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(contig, nc.range.clone(), counts);
            let reference = self.refpred.results();

            // Find loci that must be retained
            let mustloci = if let Some(ref r) = self.retainer { r.filter(contig, nc.range.clone()) } else { vec![] };

            for (strand, counts) in
                [(Strand::Forward, nc.forward), (Strand::Reverse, nc.reverse), (Strand::Unknown, nc.unstranded)]
            {
                if let Some(cnt) = counts {
                    let (stretain, stother) = self.process(
                        contig,
                        nc.range.clone(),
                        strand,
                        &mustloci,
                        cnt,
                        reference.reference,
                        reference.predicted,
                    );
                    stretain.map(|x| retained.push(x));
                    stother.map(|x| other.push(x));
                }
            }
        }

        MismatchesIntermediate { retained, other }
    }
}

struct HelperBuilder {
    pub position: Vec<Position>,
    pub refnuc: Vec<Nucleotide>,
    pub prednuc: Vec<Nucleotide>,
    pub seqnuc: Vec<NucCounts>,
}

impl HelperBuilder {
    fn new(capacity: usize) -> Self {
        Self {
            position: Vec::with_capacity(capacity),
            refnuc: Vec::with_capacity(capacity),
            prednuc: Vec::with_capacity(capacity),
            seqnuc: Vec::with_capacity(capacity),
        }
    }

    #[inline]
    fn add(&mut self, pos: Position, refnuc: Nucleotide, prednuc: Nucleotide, seqnuc: NucCounts) {
        self.position.push(pos);
        self.refnuc.push(refnuc);
        self.prednuc.push(prednuc);
        self.seqnuc.push(seqnuc);
    }

    fn finalize(self, contig: String, strand: Strand) -> Option<REATBatchedSiteMismatches> {
        if self.position.is_empty() {
            return None;
        }

        Some(REATBatchedSiteMismatches::new(contig, strand, self.position, self.refnuc, self.prednuc, self.seqnuc))
    }
}
