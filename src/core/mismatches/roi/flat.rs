use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::MismatchesSummary;
use crate::core::workload::ROI;

pub use super::ROIMismatches;

#[derive(Clone, Debug, Constructor)]
pub struct REATROIMismatches {
    roi: ROI,
    strand: Strand,
    coverage: u32,
    masked: u32,
    prednuc: NucCounts,
    mismatches: MismatchesSummary,
}

impl ROIMismatches for REATROIMismatches {
    fn roi(&self) -> &ROI {
        &self.roi
    }

    fn strand(&self) -> Strand {
        self.strand
    }

    fn coverage(&self) -> u32 {
        self.coverage
    }

    fn prednuc(&self) -> NucCounts {
        self.prednuc
    }

    fn mismatches(&self) -> &MismatchesSummary {
        &self.mismatches
    }

    fn masked(&self) -> u32 {
        self.masked
    }
}
