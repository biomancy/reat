use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::ROIRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::roi::{MismatchesSummary, REATBatchedROIMismatches, ROIMismatchesPreview};
use crate::core::mismatches::site::REATBatchedSiteMismatches;
use crate::core::mismatches::{BatchedMismatchesBuilder, MismatchesIntermediate};
use crate::core::refpred::{RefEngine, RefEngineResult};
use crate::core::rpileup::ncounters::cnt::ROINucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::workload::ROI;

pub struct REATROIMismatchesBuilder<RE: RefEngine, RR: ROIRetainer, MP: MismatchesPreFilter<ROIMismatchesPreview>> {
    buffer: Vec<NucCounts>,
    refpred: RE,
    retainer: Option<RR>,
    prefilter: Option<MP>,
}

impl<'a, RE, RR, MP> REATROIMismatchesBuilder<RE, RR, MP>
where
    RE: RefEngine,
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIMismatchesPreview>,
{
    pub fn new(refpred: RE, retainer: Option<RR>, prefilter: Option<MP>) -> Self {
        Self { buffer: vec![], refpred, retainer, prefilter }
    }

    fn process(
        &self,
        contig: &'a str,
        ncrange: &Range<Position>,
        info: &'a <ROINucCounts<'a> as AggregatedNucCounts<'a>>::ItemInfo,
        cnt: &'a [NucCounts],
        refpred: &RefEngineResult<'_>,
        retain: &mut HelperBuilder<'a>,
        other: &mut HelperBuilder<'a>,
    ) {
        let roi = info.roi;
        debug_assert!(contig == info.roi.contig());
        debug_assert!(ncrange.start <= roi.range().start && ncrange.end >= roi.range().end);

        // Get mismatches
        let (cnts, mismatches) = self.summarize(info.roi, ncrange.clone(), refpred.predicted, cnt);
        if self.retainer.as_ref().map_or(false, |x| x.filter(contig, &roi.range(), roi.name())) {
            // Must be retained
            retain.add(roi, info.coverage, cnts, mismatches);
        } else if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&mismatches)) {
            // Must be other
            other.add(roi, info.coverage, cnts, mismatches);
        }
    }

    fn summarize(
        &self,
        roi: &'a ROI,
        nucrange: Range<Position>,
        prednuc: &'a [Nucleotide],
        seqnuc: &'a [NucCounts],
    ) -> (NucCounts, MismatchesSummary) {
        let mut mismatches = MismatchesSummary::zeros();
        let mut cnts = NucCounts::zeros();

        let start = nucrange.start;
        for piece in roi.include() {
            let idx = (piece.start - start) as usize..(piece.end - start) as usize;
            cnts.increment(&prednuc[idx.clone()]);
            mismatches.increment(&prednuc[idx.clone()], &seqnuc[idx])
        }

        (cnts, mismatches)
    }

    #[inline]
    fn allocate(&self, nc: &ROINucCounts<'a>, strand: Strand) -> (HelperBuilder<'a>, HelperBuilder<'a>) {
        let iter = nc.items().iter();
        let maxitems = match strand {
            Strand::Forward => iter.map(|x| x.forward.map_or(0, |x| x.len())).sum::<usize>(),
            Strand::Reverse => iter.map(|x| x.reverse.map_or(0, |x| x.len())).sum::<usize>(),
            Strand::Unknown => iter.map(|x| x.unstranded.map_or(0, |x| x.len())).sum::<usize>(),
        };
        (HelperBuilder::new(maxitems), HelperBuilder::new(maxitems))
    }
}

impl<'a, RE, RR, MP> BatchedMismatchesBuilder<'a, ROINucCounts<'a>> for REATROIMismatchesBuilder<RE, RR, MP>
where
    RE: RefEngine,
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIMismatchesPreview>,
{
    type Mismatches = REATBatchedROIMismatches;

    fn build(&mut self, nc: ROINucCounts<'a>) -> MismatchesIntermediate<Self::Mismatches> {
        let (mut fwdretain, mut fwdother) = self.allocate(&nc, Strand::Forward);
        let (mut revretain, mut revother) = self.allocate(&nc, Strand::Reverse);
        let (mut unkretain, mut unkother) = self.allocate(&nc, Strand::Unknown);

        let (contig, items) = nc.consume();

        // 1 Item = 1 ROI
        for nc in items.iter() {
            debug_assert!(nc.forward.is_some() || nc.reverse.is_some() || nc.unstranded.is_some());

            // Predict the reference
            let counts = nc.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(contig, nc.range.clone(), counts);
            let refpred = self.refpred.results();

            // Process the counts
            nc.forward.map(|x| self.process(contig, &nc.range, &nc.info, x, &refpred, &mut fwdretain, &mut fwdother));
            nc.reverse.map(|x| self.process(contig, &nc.range, &nc.info, x, &refpred, &mut revretain, &mut revother));
            nc.unstranded
                .map(|x| self.process(contig, &nc.range, &nc.info, x, &refpred, &mut unkretain, &mut unkother));
        }

        MismatchesIntermediate {
            retained: HelperBuilder::finalize(contig, fwdretain, revretain, unkretain),
            other: HelperBuilder::finalize(contig, fwdother, revother, unkother),
        }
    }
}

struct HelperBuilder<'a> {
    pub rois: Vec<&'a ROI>,
    pub coverage: Vec<u32>,
    pub prednuc: Vec<NucCounts>,
    pub mismatches: Vec<MismatchesSummary>,
}

impl<'a> HelperBuilder<'a> {
    fn new(capacity: usize) -> Self {
        Self {
            rois: Vec::with_capacity(capacity),
            coverage: Vec::with_capacity(capacity),
            prednuc: Vec::with_capacity(capacity),
            mismatches: Vec::with_capacity(capacity),
        }
    }

    #[inline]
    fn add(&mut self, roi: &'a ROI, coverage: u32, prednuc: NucCounts, mm: MismatchesSummary) {
        self.rois.push(roi);
        self.coverage.push(coverage);
        self.prednuc.push(prednuc);
        self.mismatches.push(mm);
    }

    fn finalize(
        contig: &str,
        fwd: HelperBuilder,
        rev: HelperBuilder,
        unk: HelperBuilder,
    ) -> Vec<REATBatchedROIMismatches> {
        let mut result = Vec::with_capacity(3);

        for (strand, item) in [(Strand::Forward, fwd), (Strand::Reverse, rev), (Strand::Unknown, unk)] {
            if !item.rois.is_empty() {
                result.push(REATBatchedROIMismatches::new(
                    contig.into(),
                    strand,
                    item.rois,
                    item.coverage,
                    item.prednuc,
                    item.mismatches,
                ))
            }
        }
        result
    }
}
