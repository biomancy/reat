use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;
use itertools::{chain, izip};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::SitesRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::site::REATSiteMismatchesVec;
use crate::core::mismatches::Context;
use crate::core::refpred::{RefEngine, RefEngineResult};
use crate::core::rpileup::ncounters::cnt::IntervalNucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::strandutil::StrandedData;

use super::super::MismatchesBuilder;
use super::SiteMismatchesPreview;

#[derive(Clone)]
pub struct REATSiteMismatchesBuilder<RE, SR, MP> {
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
    pub fn new(maxsize: usize, refpred: RE, retainer: Option<SR>, prefilter: Option<MP>) -> Self {
        Self { buffer: Vec::with_capacity(maxsize), refpred, retainer, prefilter }
    }

    fn process(
        &self,
        retained: &[Range<Position>],
        cntrange: Range<Position>,
        counts: &[NucCounts],
        refngn: &RefEngineResult,
        retbuilder: &mut Helper,
        othbuilder: &mut Helper,
    ) {
        debug_assert_eq!(cntrange.end - cntrange.start, counts.len() as Position);
        debug_assert_eq!(counts.len(), refngn.reference.len());
        debug_assert_eq!(counts.len(), refngn.predicted.len());
        debug_assert!(retained.iter().all(|x| cntrange.contains(&x.start) && cntrange.contains(&x.end)));

        let retsize = retained.iter().map(|x| x.end - x.start).sum::<Position>() as usize;
        debug_assert!(retsize <= counts.len());

        // Retain iterator
        let mut reiter = retained.iter();
        let mut retrange = reiter.next();

        for (loc, (&cnt, &refn, &predn)) in izip!(counts, refngn.reference, refngn.predicted).enumerate() {
            let loc = loc as Position + cntrange.start;
            // Do we need to move the iterator?
            if !retrange.map_or(true, |x| x.end <= loc) {
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
    }
}

impl<'a, RE, SR, MP> MismatchesBuilder<'a, IntervalNucCounts<'a>> for REATSiteMismatchesBuilder<RE, SR, MP>
where
    RE: RefEngine,
    SR: SitesRetainer,
    MP: MismatchesPreFilter<SiteMismatchesPreview>,
{
    type Mismatches = REATSiteMismatchesVec;

    fn build(&mut self, nc: IntervalNucCounts<'a>) -> Context<Self::Mismatches> {
        // let capacity = nc.items().iter().map(|x| x.range.end - x.range.start).sum();
        let (mut retained, mut other) = (StrandedData::default(), StrandedData::default());

        let (contig, items) = nc.consume();
        for nc in items {
            // Predict the reference
            let counts = nc.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(contig, nc.range.clone(), counts);
            let reference = self.refpred.results();

            // Find loci that must be retained
            let mustloci = self.retainer.as_ref().map_or(vec![], |r| r.retained(contig, nc.range.clone()));
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                nc.counts[strand].map(|cnt| {
                    debug_assert!(nc.mapped[strand] > 0);
                    self.process(
                        &mustloci,
                        nc.range.clone(),
                        cnt,
                        &reference,
                        &mut retained[strand],
                        &mut other[strand],
                    );
                });
            }
        }
        todo!()
        // Context { retained: Helper::finalize(), items: Default::default() }
    }
}

#[derive(Default)]
struct Helper {
    pub position: Vec<Position>,
    pub refnuc: Vec<Nucleotide>,
    pub prednuc: Vec<Nucleotide>,
    pub seqnuc: Vec<NucCounts>,
}

impl Helper {
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

    // fn finalize(self, contig: String, strand: Strand) -> Vec<REATSiteMismatchesVec> {
    //     if self.position.is_empty() {
    //         return vec![];
    //     }
    //
    //     Some(REATSiteMismatchesVec::new(contig, strand, self.position, self.refnuc, self.prednuc, self.seqnuc))
    // }
}
