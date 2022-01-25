use bio_types::genome::Interval;
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::borrowed::RefROIMismatches;
use crate::core::mismatches::roi::{MismatchesSummary, ROIMismatches};
use crate::core::workload::ROI;

#[derive(Clone, Debug, Constructor)]
pub struct OwnedROIMismatches {
    roi: ROI,
    strand: Strand,
    masked: u32,
    coverage: u32,
    sequenced: NucCounts,
    mismatches: MismatchesSummary,
}

impl ROIMismatches for OwnedROIMismatches {
    fn roi(&self) -> &ROI {
        &self.roi
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
        &self.sequenced
    }

    fn mismatches(&self) -> &MismatchesSummary {
        &self.mismatches
    }
}

impl From<RefROIMismatches<'_>> for OwnedROIMismatches {
    fn from(x: RefROIMismatches<'_>) -> Self {
        let (roi, strand, masked, coverage, sequenced, mismatches) = x.dissolve();
        Self { roi: roi.clone(), strand, masked, coverage, sequenced, mismatches }
    }
}
