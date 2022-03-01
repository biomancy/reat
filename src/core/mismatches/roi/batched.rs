use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::flat::REATROIMismatches;
use crate::core::mismatches::roi::MismatchesSummary;
use crate::core::mismatches::{BatchedMismatches, Restranded};
use crate::core::workload::{ROIWorkload, ROI};

use super::BatchedROIMismatches;

pub struct REATBatchedROIMismatches {
    contig: String,
    range: Range<Position>,
    strand: Strand,
    names: Vec<String>,
    premasked: Vec<Range<Position>>,
    rois: Vec<Range<Position>>,
    pieces: Vec<Vec<Range<Position>>>,
    coverage: Vec<u32>,
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
        // Calculate mismatches only over the retained subintervals
        // debug_assert!(
        //     (roi.range().end - roi.range().start) as usize == prednuc.len() && prednuc.len() == sequenced.len()
        // );
        // let (mut cnts, mut mismatches) = (NucCounts::zeros(), MismatchesSummary::zeros());
        // let start = roi.range().start;
        // for piece in roi.include() {
        //     let idx = (piece.start - start) as usize..(piece.end - start) as usize;
        //     cnts.increment(&prednuc[idx.clone()]);
        //     mismatches.increment(&prednuc[idx.clone()], &sequenced[idx])
        // }

        todo!()
    }
}

impl AbstractInterval for REATBatchedROIMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

impl BatchedMismatches for REATBatchedROIMismatches {
    type Flattened = REATROIMismatches;

    fn strand(&self) -> Strand {
        self.strand
    }

    fn filter(self, mask: Vec<bool>) -> Self {
        todo!()
    }

    fn restrand(self, strands: Vec<Strand>) -> Restranded<Self> {
        todo!()
    }

    fn restrand_all(&mut self, strand: Strand) {
        self.strand = strand;
    }

    fn flatten(self) -> Vec<Self::Flattened> {
        todo!()
    }
}

impl BatchedROIMismatches for REATBatchedROIMismatches {
    fn rois(&self) -> &[Range<Position>] {
        &self.rois
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
