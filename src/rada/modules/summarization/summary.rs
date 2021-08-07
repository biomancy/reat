use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;
use derive_more::Constructor;

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::ReqNucleotide;

#[derive(Constructor)]
#[allow(non_snake_case)]
pub struct MismatchSummary {
    pub A: LocusCounts,
    pub C: LocusCounts,
    pub G: LocusCounts,
    pub T: LocusCounts,
}

impl MismatchSummary {
    pub fn zeros() -> Self {
        MismatchSummary {
            A: LocusCounts::zeros(),
            C: LocusCounts::zeros(),
            G: LocusCounts::zeros(),
            T: LocusCounts::zeros(),
        }
    }

    pub fn coverage(&self) -> u32 {
        self.A.coverage() + self.C.coverage() + self.G.coverage() + self.T.coverage()
    }

    pub fn matches(&self) -> u32 {
        self.A.A + self.C.C + self.G.G + self.T.T
    }

    pub fn mismatches(&self) -> u32 {
        self.A.mismatches(&ReqNucleotide::A)
            + self.C.mismatches(&ReqNucleotide::C)
            + self.G.mismatches(&ReqNucleotide::G)
            + self.T.mismatches(&ReqNucleotide::T)
    }
}

pub struct IntervalSummary {
    pub interval: Interval,
    pub strand: Strand,
    pub name: String,
    pub mismatches: MismatchSummary,
}

impl IntervalSummary {
    pub fn zeros(interval: Interval, strand: Strand, name: String) -> Self {
        IntervalSummary {
            interval,
            strand,
            name,
            mismatches: MismatchSummary::zeros(),
        }
    }
}

pub struct LocusSummary {
    pub locus: Locus,
    pub strand: Strand,
    pub refnuc: ReqNucleotide,
    pub sequenced: LocusCounts,
}

impl LocusSummary {
    pub fn zeros(locus: Locus, strand: Strand, refnuc: ReqNucleotide) -> Self {
        LocusSummary {
            locus,
            strand,
            refnuc,
            sequenced: LocusCounts::zeros(),
        }
    }
}
