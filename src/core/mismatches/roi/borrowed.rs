use bio_types::genome::Interval;
use bio_types::strand::{ReqStrand, Strand};
use derive_getters::Dissolve;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{MismatchesSummary, OwnedROIMismatches, ROIMismatches};
use crate::core::mismatches::IntermediateMismatches;

#[derive(Clone, Debug, Dissolve)]
pub struct BorrowedROIMismatches<'a> {
    interval: &'a Interval,
    strand: Strand,
    name: &'a String,
    coverage: u32,
    sequence: NucCounts,
    mismatches: MismatchesSummary,
}

impl ROIMismatches for BorrowedROIMismatches<'_> {
    fn interval(&self) -> &Interval {
        &self.interval
    }

    fn strand(&self) -> &Strand {
        &self.strand
    }

    fn name(&self) -> &String {
        &self.name
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

impl<'a> BorrowedROIMismatches<'a> {
    pub fn new(
        interval: &'a Interval,
        strand: Strand,
        name: &'a String,
        coverage: u32,
        reference: &'a [Nucleotide],
        sequenced: &'a [NucCounts],
    ) -> BorrowedROIMismatches<'a> {
        Self {
            interval,
            strand,
            name,
            coverage,
            sequence: NucCounts::from_sequence(reference),
            mismatches: MismatchesSummary::from_ncounts(reference, sequenced),
        }
    }
}

impl IntermediateMismatches for BorrowedROIMismatches<'_> {
    type Finalized = OwnedROIMismatches;

    fn finalize(self) -> Self::Finalized {
        self.into()
    }

    fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }
}
