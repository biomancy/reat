use derive_more::Constructor;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::Nucleotide;

use super::RefNucPredictor;
use derive_getters::Getters;

#[derive(Constructor, Getters, Copy, Clone)]
pub struct RefNucPredByHeurisitc {
    mincoverage: u32,
    freqthr: f32,
}

impl RefNucPredictor for RefNucPredByHeurisitc {
    fn predict(&self, assembly: &Nucleotide, sequenced: &LocusCounts) -> Nucleotide {
        let coverage = sequenced.coverage();

        if coverage >= self.mincoverage {
            let (nuc, counts) = sequenced.mostfreq();
            if *counts as f32 / coverage as f32 >= self.freqthr {
                return nuc.into();
            }
        }
        *assembly
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn predict() {
        let sequenced = LocusCounts { A: 1, C: 2, G: 3, T: 4 };
        let assembly = Nucleotide::A;
        for (mincoverage, freqthr, result) in [
            (100, 0.0, assembly),
            (0, 1.0, assembly),
            (0, 0.41, assembly),
            (0, 0.0, Nucleotide::T),
            (0, 0.4, Nucleotide::T),
            (4, 0.4, Nucleotide::T),
        ] {
            let dummy = RefNucPredByHeurisitc { mincoverage, freqthr };
            assert_eq!(dummy.predict(&assembly, &sequenced), result);
        }
    }
}
