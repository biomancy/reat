

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::prefilters::retain::ROIRetainer;
use crate::core::mismatches::prefilters::MismatchesPreFilter;
use crate::core::mismatches::roi::{NucMismatches, ROIData, ROIDataVec, ROIMismatchesVec};
use crate::core::mismatches::{Batch, Builder};
use crate::core::refpred::{RefEngine, RefEngineResult};
use crate::core::rpileup::ncounter::NucCounterResult;
use crate::core::strandutil::Stranded;
use crate::core::workload::ROI;

#[derive(Clone)]
pub struct ROIMismatchesBuilder<RE: RefEngine, RR: ROIRetainer, MP: MismatchesPreFilter<ROIData>> {
    buffer: Vec<NucCounts>,
    refpred: RE,
    retainer: Option<RR>,
    prefilter: Option<MP>,
}

impl<'a, RE, RR, MP> ROIMismatchesBuilder<RE, RR, MP>
where
    RE: RefEngine,
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIData>,
{
    pub fn new(maxsize: usize, refpred: RE, retainer: Option<RR>, prefilter: Option<MP>) -> Self {
        Self { buffer: Vec::with_capacity(maxsize), refpred, retainer, prefilter }
    }

    fn process(
        &self,
        cntstart: Position,
        cnts: &'a [NucCounts],
        refpred: &RefEngineResult<'_>,
        roi: &'a ROI,
        coverage: u32,
        retain: &mut ROIDataVec,
        other: &mut ROIDataVec,
    ) {
        // Get mismatches
        let (prednuc, mismatches) = self.summarize(roi, cntstart, refpred.predicted, cnts);
        let record = ROIData { roi: roi.into(), coverage, prednuc, mismatches };
        if self.retainer.as_ref().map_or(false, |x| x.retained(roi.contig(), &roi.range(), roi.strand(), roi.name())) {
            // Must be retained
            retain.push(record);
        } else if self.prefilter.as_ref().map_or(true, |x| x.is_ok(&record)) {
            // Must be other
            other.push(record);
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

    fn finalize(contig: &str, item: Stranded<ROIDataVec>) -> Stranded<ROIMismatchesVec> {
        item.into(|x, strand| ROIMismatchesVec::new(contig.to_owned(), strand, x))
    }

    #[inline]
    fn size_hint(&self, nc: &NucCounterResult<'a, &'a ROI>) -> Stranded<usize> {
        let mut size: Stranded<usize> = Stranded::default();
        for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
            size[strand] = nc.cnts.iter().map(|item| item.cnts[strand].is_some() as usize).sum()
        }
        size
    }
}

impl<'a, RE, RR, MP> Builder<'a> for ROIMismatchesBuilder<RE, RR, MP>
where
    RE: RefEngine,
    RR: ROIRetainer,
    MP: MismatchesPreFilter<ROIData>,
{
    type Out = ROIMismatchesVec;
    type SourceCounts = NucCounterResult<'a, &'a ROI>;

    fn build(&mut self, nc: Self::SourceCounts) -> Batch<Self::Out> {
        let contig = nc.contig.to_owned();

        // Pre-allocate results
        let _hint = self.size_hint(&nc);
        let mut items = Stranded::with_fn(|strnd| {
            ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::new())
            // ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::with_capacity(hint[strnd]))
        });
        let mut retained = Stranded::with_fn(|strnd| {
            ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::new())
            // ROIMismatchesVec::new(contig.clone(), strnd, ROIDataVec::with_capacity(hint[strnd] / 10))
        });

        for item in nc.cnts.into_iter() {
            // debug_assert!(item.coverage.forward + item.coverage.reverse + item.coverage.unknown > 0);

            // Predict the reference
            let counts = item.seqnuc(&mut self.buffer).unwrap_or(&self.buffer);
            self.refpred.run(&contig, item.range.clone(), counts);
            let refpred = self.refpred.results();

            // Process the counts
            for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
                if let Some(cnts) = item.cnts[strand] {
                    self.process(
                        item.range.start,
                        cnts,
                        &refpred,
                        item.data,
                        item.coverage[strand],
                        &mut retained[strand].data,
                        &mut items[strand].data,
                    );
                }
            }
        }
        Batch { contig, mapped: nc.mapped, retained, items }
    }
}
