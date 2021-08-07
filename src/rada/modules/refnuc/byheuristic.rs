

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::{Nucleotide, ReqNucleotide};

use super::RefNucPredictor;

pub struct RefNucPredByHeurisitc {
    mincoverage: u32,
    freqthr: f32,
}

impl RefNucPredByHeurisitc {
    pub fn new(mincoverage: u32, freqthr: f32) -> Self {
        RefNucPredByHeurisitc {
            mincoverage,
            freqthr,
        }
    }
}

impl RefNucPredictor for RefNucPredByHeurisitc {
    fn predict(
        &self,
        assembly: Nucleotide,
        _: &[ReqNucleotide],
        sequenced: &LocusCounts,
    ) -> ReqNucleotide {
        let coverage = sequenced.coverage();

        if coverage > self.mincoverage {
            let coverage = coverage as f32;
            if (sequenced.A as f32 / coverage) > self.freqthr {
                return ReqNucleotide::A;
            } else if (sequenced.T as f32 / coverage) > self.freqthr {
                return ReqNucleotide::T;
            } else if (sequenced.G as f32 / coverage) > self.freqthr {
                return ReqNucleotide::G;
            } else if (sequenced.C as f32 / coverage) > self.freqthr {
                return ReqNucleotide::C;
            }
        }

        if assembly == Nucleotide::Unknown {
            sequenced.mostfreq().0
        } else {
            assembly.into()
        }
    }
}
