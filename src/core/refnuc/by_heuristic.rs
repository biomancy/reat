use derive_getters::Getters;
use derive_more::Constructor;

use crate::core::counting::LocusCounts;
use crate::core::dna::{Nucleotide, ReqNucleotide};

use super::RefNucPredictor;

#[derive(Constructor, Getters, Copy, Clone)]
pub struct RefNucPredByHeurisitc {
    mincoverage: u32,
    freqthr: f32,
    skip_hyperediting: bool,
}

impl RefNucPredictor for RefNucPredByHeurisitc {
    fn predict(&self, assembly: &Nucleotide, sequenced: &LocusCounts) -> Nucleotide {
        let coverage = sequenced.coverage();

        // if coverage is sufficient
        if coverage >= self.mincoverage {
            let (nuc, counts) = sequenced.mostfreq();
            // and the most abundant nucleotide is frequent enough
            if *counts as f32 / coverage as f32 >= self.freqthr {
                // check for a potential hyper-editing site
                let skip_hyperediting = self.skip_hyperediting
                    && ((*assembly == Nucleotide::A && nuc == ReqNucleotide::G)
                        || (*assembly == Nucleotide::T && nuc == ReqNucleotide::C));

                if skip_hyperediting {
                    return *assembly;
                } else {
                    return nuc.into();
                };
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
            let dummy = RefNucPredByHeurisitc { mincoverage, freqthr, skip_hyperediting: false };
            assert_eq!(dummy.predict(&assembly, &sequenced), result);
        }
    }

    #[test]
    fn skip_hyper_editing() {
        let sequenced = LocusCounts { A: 1, C: 0, G: 99, T: 0 };
        for (expect, toskip) in [(Nucleotide::A, true), (Nucleotide::G, false)] {
            let dummy = RefNucPredByHeurisitc { mincoverage: 0, freqthr: 0f32, skip_hyperediting: toskip };
            assert_eq!(dummy.predict(&Nucleotide::A, &sequenced), expect);
        }

        let sequenced = LocusCounts { A: 0, C: 3, G: 0, T: 1 };
        for (expect, toskip) in [(Nucleotide::T, true), (Nucleotide::C, false)] {
            let dummy = RefNucPredByHeurisitc { mincoverage: 0, freqthr: 0f32, skip_hyperediting: toskip };
            assert_eq!(dummy.predict(&Nucleotide::T, &sequenced), expect);
        }
    }
}
