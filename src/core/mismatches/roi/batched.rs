use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use bio_types::strand::Strand;
use itertools::{izip, Itertools};

use crate::core::dna::NucCounts;
use crate::core::mismatches::roi::flat::REATROIMismatches;
use crate::core::mismatches::roi::MismatchesSummary;
use crate::core::mismatches::utils::{maskvec, select_strands};
use crate::core::mismatches::{BatchedMismatches, StrandingCounts};
use crate::core::strandutil::StrandedData;
use crate::core::workload::ROI;

use super::BatchedROIMismatches;

pub struct REATBatchedROIMismatches {
    contig: String,
    strand: Strand,
    names: Vec<String>,
    annostrand: Vec<Strand>,
    coverage: Vec<u32>,
    premasked: Vec<Range<Position>>,
    postmasked: Vec<Range<Position>>,
    subintervals: Vec<Vec<Range<Position>>>,
    prednuc: Vec<NucCounts>,
    mismatches: Vec<MismatchesSummary>,
}

impl REATBatchedROIMismatches {
    pub fn new(
        contig: String,
        strand: Strand,
        rois: Vec<&ROI>,
        coverage: Vec<u32>,
        prednuc: Vec<NucCounts>,
        mismatches: Vec<MismatchesSummary>,
    ) -> Self {
        debug_assert!(!rois.is_empty());
        debug_assert!([rois.len(), coverage.len(), prednuc.len(), mismatches.len()].iter().all_equal());
        // ROIs must be ordered
        debug_assert!(rois.windows(2).all(|x| x[0].range().start <= x[1].range().start));
        let items = rois.len();
        let (
            mut names,
            mut strands,
            mut premasked,
            mut postmasked,
            mut subintervals,
            mut roicov,
            mut roiprednuc,
            mut roimm,
        ) = (
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
            Vec::with_capacity(items),
        );

        for (roi, coverage, prednuc, mm) in izip!(rois, coverage, prednuc, mismatches) {
            names.push(roi.name().clone());
            strands.push(*roi.strand());
            premasked.push(roi.interval().range());
            postmasked.push(roi.range());
            subintervals.push(roi.subintervals().clone());
            roicov.push(coverage);
            roiprednuc.push(prednuc);
            roimm.push(mm)
        }

        Self {
            contig,
            strand,
            names,
            annostrand: strands,
            coverage: roicov,
            premasked,
            postmasked,
            subintervals,
            prednuc: roiprednuc,
            mismatches: roimm,
        }
    }
}

impl AbstractInterval for REATBatchedROIMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.postmasked.first().unwrap().start..self.postmasked.last().unwrap().end
    }
}

impl BatchedMismatches for REATBatchedROIMismatches {
    type Flattened = REATROIMismatches;

    fn trstrand(&self) -> Strand {
        self.strand
    }

    fn filter(mut self, mask: Vec<bool>) -> Self {
        self.names = maskvec(self.names, &mask);
        self.annostrand = maskvec(self.annostrand, &mask);
        self.coverage = maskvec(self.coverage, &mask);
        self.premasked = maskvec(self.premasked, &mask);
        self.postmasked = maskvec(self.postmasked, &mask);
        self.subintervals = maskvec(self.subintervals, &mask);
        self.prednuc = maskvec(self.prednuc, &mask);
        self.mismatches = maskvec(self.mismatches, &mask);
        self
    }

    fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>> {
        let names = select_strands(self.names, &strands, &cnts);
        let annostrand = select_strands(self.annostrand, &strands, &cnts);
        let coverage = select_strands(self.coverage, &strands, &cnts);
        let premasked = select_strands(self.premasked, &strands, &cnts);
        let postmasked = select_strands(self.postmasked, &strands, &cnts);
        let subintervals = select_strands(self.subintervals, &strands, &cnts);
        let prednuc = select_strands(self.prednuc, &strands, &cnts);
        let mismatches = select_strands(self.mismatches, &strands, &cnts);

        let mut result = StrandedData { unknown: None, forward: None, reverse: None };
        if cnts.forward > 0 {
            result.forward = Some(Self {
                contig: self.contig.to_owned(),
                strand: Strand::Forward,
                names: names.0,
                annostrand: annostrand.0,
                coverage: coverage.0,
                premasked: premasked.0,
                postmasked: postmasked.0,
                subintervals: subintervals.0,
                prednuc: prednuc.0,
                mismatches: mismatches.0,
            })
        }
        if cnts.reverse > 0 {
            result.reverse = Some(Self {
                contig: self.contig.to_owned(),
                strand: Strand::Reverse,
                names: names.1,
                annostrand: annostrand.1,
                coverage: coverage.1,
                premasked: premasked.1,
                postmasked: postmasked.1,
                subintervals: subintervals.1,
                prednuc: prednuc.1,
                mismatches: mismatches.1,
            })
        }
        if cnts.unknown > 0 {
            result.unknown = Some(Self {
                contig: self.contig.to_owned(),
                strand: Strand::Unknown,
                names: names.2,
                annostrand: annostrand.2,
                coverage: coverage.2,
                premasked: premasked.2,
                postmasked: postmasked.2,
                subintervals: subintervals.2,
                prednuc: prednuc.2,
                mismatches: mismatches.2,
            })
        }
        result
    }

    fn restrand_all(&mut self, strand: Strand) {
        self.strand = strand;
    }

    fn flatten(self) -> Vec<Self::Flattened> {
        let mut result = Vec::with_capacity(self.names.len());
        for (name, annostrand, coverage, premasked, subintervals, prednuc, mismatches) in izip!(
            self.names,
            self.annostrand,
            self.coverage,
            self.premasked,
            self.subintervals,
            self.prednuc,
            self.mismatches
        ) {
            let roi = Interval::new(self.contig.to_owned(), premasked);
            let roi = ROI::new(roi, name, annostrand, subintervals);
            result.push(REATROIMismatches::new(roi, self.strand, coverage, prednuc, mismatches))
        }
        result
    }
}

impl BatchedROIMismatches for REATBatchedROIMismatches {
    fn premasked(&self) -> &[Range<Position>] {
        &self.premasked
    }

    fn postmasked(&self) -> &[Range<Position>] {
        &self.postmasked
    }

    fn subintervals(&self) -> &[Vec<Range<Position>>] {
        &self.subintervals
    }

    fn coverage(&self) -> &[u32] {
        &self.coverage
    }

    fn prednuc(&self) -> &[NucCounts] {
        &self.prednuc
    }

    fn mismatches(&self) -> &[MismatchesSummary] {
        &self.mismatches
    }
}
