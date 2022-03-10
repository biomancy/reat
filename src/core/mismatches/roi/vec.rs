use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::Strand;
use itertools::{izip, Itertools};

use crate::core::dna::NucCounts;
use crate::core::mismatches::roi::flat::REATROIMismatches;
use crate::core::mismatches::roi::NucMismatches;
use crate::core::mismatches::utils::{maskvec, select_strands};
use crate::core::mismatches::{MismatchesVec, StrandingCounts};
use crate::core::strandutil::StrandedData;
use crate::core::workload::ROI;

use super::ROIMismatchesVec;

pub struct REATROIMismatchesVec {
    contig: String,
    trstrand: Strand,
    // ROI data
    roi_names: Vec<String>,
    roi_strands: Vec<Strand>,
    roi_premasked: Vec<Range<Position>>,
    roi_postmasked: Vec<Range<Position>>,
    roi_subintervals: Vec<Vec<Range<Position>>>,
    // Calculated stats
    prednuc: Vec<NucCounts>,
    coverage: Vec<u32>,
    mismatches: Vec<NucMismatches>,
}

impl REATROIMismatchesVec {
    pub fn new(
        contig: String,
        trstrand: Strand,
        rois: Vec<&ROI>,
        coverage: Vec<u32>,
        prednuc: Vec<NucCounts>,
        mismatches: Vec<NucMismatches>,
    ) -> Self {
        debug_assert!(!rois.is_empty());
        debug_assert!([rois.len(), coverage.len(), prednuc.len(), mismatches.len()].iter().all_equal());
        // ROIs must be ordered
        debug_assert!(rois.windows(2).all(|x| x[0].range().start <= x[1].range().start));

        let items = rois.len();
        let mut results = Self::with_capacity(contig, trstrand, items);

        for (roi, coverage, prednuc, mm) in izip!(rois, coverage, prednuc, mismatches) {
            debug_assert_eq!(roi.contig(), results.contig);
            results.mismatches.push(mm);
            results.coverage.push(coverage);
            results.prednuc.push(prednuc);

            // ROI data
            results.roi_premasked.push(roi.premasked());
            results.roi_strands.push(roi.strand());
            results.roi_names.push(roi.name().into());
            results.roi_postmasked.push(roi.postmasked());
            results.roi_subintervals.push(roi.subintervals().to_vec());
        }
        results
    }

    fn with_capacity(contig: String, trstrand: Strand, capacity: usize) -> Self {
        Self {
            contig,
            trstrand,
            roi_names: Vec::with_capacity(capacity),
            roi_strands: Vec::with_capacity(capacity),
            roi_premasked: Vec::with_capacity(capacity),
            roi_postmasked: Vec::with_capacity(capacity),
            roi_subintervals: Vec::with_capacity(capacity),
            prednuc: Vec::with_capacity(capacity),
            coverage: Vec::with_capacity(capacity),
            mismatches: Vec::with_capacity(capacity),
        }
    }
}

impl AbstractInterval for REATROIMismatchesVec {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.roi_postmasked.first().unwrap().start..self.roi_postmasked.last().unwrap().end
    }
}

impl MismatchesVec for REATROIMismatchesVec {
    type Flat = REATROIMismatches;

    fn trstrand(&self) -> Strand {
        self.trstrand
    }

    // fn filter(mut self, mask: Vec<bool>) -> Self {
    //     self.roi_premasked = maskvec(self.roi_premasked, &mask);
    //     self.roi_strands = maskvec(self.roi_strands, &mask);
    //     self.roi_names = maskvec(self.roi_names, &mask);
    //     self.roi_subintervals = maskvec(self.roi_subintervals, &mask);
    //     self.roi_postmasked = maskvec(self.roi_postmasked, &mask);
    //
    //     self.prednuc = maskvec(self.prednuc, &mask);
    //     self.coverage = maskvec(self.coverage, &mask);
    //     self.mismatches = maskvec(self.mismatches, &mask);
    //
    //     self
    // }
    //
    // fn extend(&mut self, other: Self) {
    //     todo!()
    // }
    //
    // fn restrand(self, strands: Vec<Strand>, fwdto: &mut Self, revto: &mut Self) {
    //     todo!()
    // }

    // fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>> {
    //     let mut names = select_strands(self.roi_names, &strands, &cnts);
    //     let mut annostrand = select_strands(self.roi_strands, &strands, &cnts);
    //     let mut coverage = select_strands(self.coverage, &strands, &cnts);
    //     let mut premasked = select_strands(self.roi_premasked, &strands, &cnts);
    //     let mut postmasked = select_strands(self.roi_postmasked, &strands, &cnts);
    //     let mut subintervals = select_strands(self.roi_subintervals, &strands, &cnts);
    //     let mut prednuc = select_strands(self.prednuc, &strands, &cnts);
    //     let mut mismatches = select_strands(self.mismatches, &strands, &cnts);
    //
    //     let mut result = StrandedData { unknown: None, forward: None, reverse: None };
    //     for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
    //         if cnts[strand] > 0 {
    //             result[strand] = Some(Self {
    //                 contig: self.contig.clone(),
    //                 trstrand: strand,
    //                 roi_premasked: std::mem::take(&mut premasked[strand]),
    //                 roi_strands: std::mem::take(&mut annostrand[strand]),
    //                 roi_names: std::mem::take(&mut names[strand]),
    //                 roi_subintervals: std::mem::take(&mut subintervals[strand]),
    //                 roi_postmasked: std::mem::take(&mut postmasked[strand]),
    //                 prednuc: std::mem::take(&mut prednuc[strand]),
    //                 coverage: std::mem::take(&mut coverage[strand]),
    //                 mismatches: std::mem::take(&mut mismatches[strand]),
    //             })
    //         }
    //     }
    //     result
    // }

    // fn restrand_all(&mut self, strand: Strand) {
    //     self.trstrand = strand;
    // }

    fn flatten(self) -> Vec<Self::Flat> {
        let mut result = Vec::with_capacity(self.roi_names.len());
        for (name, annostrand, coverage, premasked, subintervals, prednuc, mismatches) in izip!(
            self.roi_names,
            self.roi_strands,
            self.coverage,
            self.roi_premasked,
            self.roi_subintervals,
            self.prednuc,
            self.mismatches
        ) {
            let roi = ROI::new(self.contig.to_owned(), premasked, subintervals, name, annostrand);
            result.push(REATROIMismatches::new(roi, self.trstrand, coverage, prednuc, mismatches))
        }
        result
    }
}

impl ROIMismatchesVec for REATROIMismatchesVec {
    fn premasked(&self) -> &[Range<Position>] {
        &self.roi_premasked
    }

    fn postmasked(&self) -> &[Range<Position>] {
        &self.roi_postmasked
    }

    fn subintervals(&self) -> &[Vec<Range<Position>>] {
        &self.roi_subintervals
    }

    fn coverage(&self) -> &[u32] {
        &self.coverage
    }

    fn prednuc(&self) -> &[NucCounts] {
        &self.prednuc
    }

    fn mismatches(&self) -> &[NucMismatches] {
        &self.mismatches
    }
}
