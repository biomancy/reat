use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::ROIRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::roi::{NucMismatches, REATROIMismatchesVec, ROIMismatchesPreview};
use crate::core::mismatches::{Context, MismatchesBuilder};
use crate::core::refpred::{RefEngine, RefEngineResult};
use crate::core::rpileup::ncounters::cnt::ROINucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCounts;
use crate::core::strandutil::StrandedData;
use crate::core::workload::ROI;

#[derive(Clone)]
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
    pub fn new(maxsize: usize, refpred: RE, retainer: Option<RR>, prefilter: Option<MP>) -> Self {
        Self { buffer: Vec::with_capacity(maxsize), refpred, retainer, prefilter }
    }

    fn process(
        &self,
        contig: &'a str,
        cntstart: Position,
        cnts: &'a [NucCounts],
        refpred: &RefEngineResult<'_>,
        roi: &'a ROI,
        coverage: u32,
        retain: &mut HelperBuilder<'a>,
        other: &mut HelperBuilder<'a>,
    ) {
        debug_assert!(contig == roi.contig());

        // Get mismatches
        let (cnts, mismatches) = self.summarize(roi, cntstart, refpred.predicted, cnts);
        if self.retainer.as_ref().map_or(false, |x| x.retained(contig, &roi.range(), roi.strand(), roi.name())) {
            // Must be retained
            retain.add(roi, coverage, cnts, mismatches);
        } else if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&mismatches)) {
            // Must be other
            other.add(roi, coverage, cnts, mismatches);
        }
    }

    fn summarize(
        &self,
        roi: &'a ROI,
        cntstart: Position,
        prednuc: &'a [Nucleotide],
        cnts: &'a [NucCounts],
    ) -> (NucCounts, NucMismatches) {
        debug_assert!(roi.range().start >= cntstart && roi.range().end <= (cntstart + cnts.len() as u64));
        let mut mismatches = NucMismatches::zeros();
        let mut nuccnts = NucCounts::zeros();

        for sub in roi.subintervals() {
            let idx = (sub.start - cntstart) as usize..(sub.end - cntstart) as usize;
            nuccnts.increment(&prednuc[idx.clone()]);
            mismatches.increment(&prednuc[idx.clone()], &cnts[idx])
        }

        (nuccnts, mismatches)
    }

    #[inline]
    fn allocate(&self, nc: &ROINucCounts<'a>) -> (StrandedData<HelperBuilder<'a>>, StrandedData<HelperBuilder<'a>>) {
        let mut allocate: StrandedData<usize> = StrandedData::default();
        for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
            allocate[strand] = nc.items().iter().map(|item| item.counts[strand].map_or(0, |cnts| cnts.len())).sum()
        }

        let helper = StrandedData {
            forward: HelperBuilder::with_capacity(allocate.forward),
            reverse: HelperBuilder::with_capacity(allocate.reverse),
            unknown: HelperBuilder::with_capacity(allocate.unknown),
        };
        (helper.clone(), helper)
    }
}

impl<'a, RE, RR, MP> MismatchesBuilder<'a, ROINucCounts<'a>> for REATROIMismatchesBuilder<RE, RR, MP>
where
    RE: RefEngine,
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIMismatchesPreview>,
{
    type Mismatches = REATROIMismatchesVec;

    fn build(&mut self, nc: ROINucCounts<'a>) -> Context<Self::Mismatches> {
        let (mut retained, mut other) = self.allocate(&nc);

        let (contig, items) = nc.consume();
        for nc in items.into_iter() {
            debug_assert!(nc.mapped.forward + nc.mapped.reverse + nc.mapped.unknown > 0);

            // Predict the reference
            let counts = nc.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(contig, nc.range.clone(), counts);
            let refpred = self.refpred.results();

            // Process the counts
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                nc.counts[strand].map(|cnts| {
                    self.process(
                        contig,
                        nc.range.start,
                        cnts,
                        &refpred,
                        nc.info.roi,
                        nc.info.coverage,
                        &mut retained[strand],
                        &mut other[strand],
                    );
                });
            }
        }
        todo!()
    }
}

#[derive(Clone)]
struct HelperBuilder<'a> {
    pub rois: Vec<&'a ROI>,
    pub coverage: Vec<u32>,
    pub prednuc: Vec<NucCounts>,
    pub mismatches: Vec<NucMismatches>,
}

impl<'a> HelperBuilder<'a> {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            rois: Vec::with_capacity(capacity),
            coverage: Vec::with_capacity(capacity),
            prednuc: Vec::with_capacity(capacity),
            mismatches: Vec::with_capacity(capacity),
        }
    }

    #[inline]
    fn add(&mut self, roi: &'a ROI, coverage: u32, prednuc: NucCounts, mm: NucMismatches) {
        self.rois.push(roi);
        self.coverage.push(coverage);
        self.prednuc.push(prednuc);
        self.mismatches.push(mm);
    }

    fn finalize(contig: &str, fwd: HelperBuilder, rev: HelperBuilder, unk: HelperBuilder) -> Vec<REATROIMismatchesVec> {
        let mut result = Vec::with_capacity(3);

        for (strand, item) in [(Strand::Forward, fwd), (Strand::Reverse, rev), (Strand::Unknown, unk)] {
            if !item.rois.is_empty() {
                result.push(REATROIMismatchesVec::new(
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
