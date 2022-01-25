use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{ReqStrand, Strand};
use derive_getters::Dissolve;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{MismatchesSummary, OwnedROIMismatches, ROIMismatches};
use crate::core::mismatches::IntermediateMismatches;
use crate::core::workload::ROI;

#[derive(Clone, Debug, Dissolve)]
pub struct RefROIMismatches<'a> {
    roi: &'a ROI,
    strand: Strand,
    masked: u32,
    coverage: u32,
    sequence: NucCounts,
    mismatches: MismatchesSummary,
}

impl ROIMismatches for RefROIMismatches<'_> {
    fn roi(&self) -> &ROI {
        self.roi
    }

    fn strand(&self) -> &Strand {
        &self.strand
    }

    fn masked(&self) -> u32 {
        self.masked
    }

    fn coverage(&self) -> u32 {
        self.coverage
    }

    fn sequence(&self) -> &NucCounts {
        &self.sequence
    }

    fn mismatches(&self) -> &MismatchesSummary {
        &self.mismatches
    }
}

impl<'a> RefROIMismatches<'a> {
    pub fn new(
        roi: &'a ROI,
        strand: Strand,
        coverage: u32,
        reference: &'a [Nucleotide],
        sequenced: &'a [NucCounts],
    ) -> RefROIMismatches<'a> {
        debug_assert!(
            (roi.range().end - roi.range().start) as usize == reference.len() && reference.len() == sequenced.len()
        );
        // Calculate mismatches only over the retained subintervals
        let (mut cnts, mut mismatches) = (NucCounts::zeros(), MismatchesSummary::zeros());
        let start = roi.range().start;
        let mut nucin = 0;
        for piece in roi.include() {
            nucin += piece.end - piece.start;

            let idx = (piece.start - start) as usize..(piece.end - start) as usize;
            cnts.increment(&reference[idx.clone()]);
            mismatches.increment(&reference[idx.clone()], &sequenced[idx])
        }

        let total = roi.original().range().end - roi.original().range().start;
        Self { roi, strand, masked: (total - nucin) as u32, coverage, sequence: cnts, mismatches }
    }
}

impl IntermediateMismatches for RefROIMismatches<'_> {
    type Finalized = OwnedROIMismatches;

    fn finalize(self) -> Self::Finalized {
        self.into()
    }

    fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }
}
