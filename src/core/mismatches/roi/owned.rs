use bio_types::genome::Interval;
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::borrowed::RefROIMismatches;
use crate::core::mismatches::roi::{MismatchesSummary, ROIMismatches};

#[derive(Clone, Debug, Constructor)]
pub struct OwnedROIMismatches {
    interval: Interval,
    strand: Strand,
    name: String,
    coverage: u32,
    sequenced: NucCounts,
    mismatches: MismatchesSummary,
}

impl ROIMismatches for OwnedROIMismatches {
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
        &self.sequenced
    }

    fn mismatches(&self) -> &MismatchesSummary {
        &self.mismatches
    }
}

impl From<RefROIMismatches<'_>> for OwnedROIMismatches {
    fn from(x: RefROIMismatches<'_>) -> Self {
        let (interval, strand, name, coverage, sequenced, mismatches) = x.dissolve();
        Self { interval: interval.clone(), strand, name: name.clone(), coverage, sequenced, mismatches }
    }
}
