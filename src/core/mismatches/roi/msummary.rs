use derive_more::{Add, AddAssign, Constructor};
use itertools::zip;

pub use crate::core::dna::ncounts::*;
use crate::core::dna::Nucleotide;

// Simple struct-level attribute is not working for some reason
mod inner {
    #![allow(non_snake_case)]
    use super::*;
    // These are counts with respect to the FORWARD strand
    #[derive(Constructor, Copy, Clone, Eq, PartialEq, Debug, Add, AddAssign)]
    pub struct NucMismatches {
        pub A: NucCounts,
        pub C: NucCounts,
        pub G: NucCounts,
        pub T: NucCounts,
    }
}
pub use inner::NucMismatches;

impl NucMismatches {
    #[inline]
    pub fn zeros() -> Self {
        NucMismatches { A: NucCounts::zeros(), C: NucCounts::zeros(), G: NucCounts::zeros(), T: NucCounts::zeros() }
    }

    pub fn increment(&mut self, sequence: &[Nucleotide], ncounts: &[NucCounts]) {
        for (seq, count) in zip(sequence, ncounts) {
            match seq {
                Nucleotide::A => {
                    self.A += *count;
                }
                Nucleotide::C => {
                    self.C += *count;
                }
                Nucleotide::G => {
                    self.G += *count;
                }
                Nucleotide::T => {
                    self.T += *count;
                }
                _ => {}
            }
        }
    }

    #[inline]
    pub fn coverage(&self) -> u32 {
        self.A.coverage() + self.C.coverage() + self.G.coverage() + self.T.coverage()
    }

    #[inline]
    pub fn mismatches(&self) -> u32 {
        self.A.mismatches(Nucleotide::A)
            + self.C.mismatches(Nucleotide::C)
            + self.G.mismatches(Nucleotide::G)
            + self.T.mismatches(Nucleotide::T)
    }

    #[inline]
    pub fn complementary(&self) -> Self {
        NucMismatches {
            A: self.T.complementary(),
            C: self.G.complementary(),
            G: self.C.complementary(),
            T: self.A.complementary(),
        }
    }
}

impl Default for NucMismatches {
    fn default() -> Self {
        NucMismatches::zeros()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fillall(value: u32, counts: &mut NucCounts) {
        counts.A = value;
        counts.C = value;
        counts.G = value;
        counts.T = value;
    }

    #[test]
    fn coverage() {
        let mut dummy: NucMismatches = Default::default();
        for num in [0, 25] {
            fillall(num, &mut dummy.A);
            fillall(num, &mut dummy.C);
            fillall(num, &mut dummy.G);
            fillall(num, &mut dummy.T);
            assert_eq!(dummy.coverage(), num * 4 * 4);
        }
    }

    #[test]
    fn match_mismatch() {
        let mut dummy: NucMismatches = Default::default();
        for (mismatched, matched) in [(0, 25), (12, 0), (1, 2)] {
            fillall(mismatched, &mut dummy.A);
            dummy.A.A = matched;
            fillall(mismatched, &mut dummy.C);
            dummy.C.C = matched;
            fillall(mismatched, &mut dummy.G);
            dummy.G.G = matched;
            fillall(mismatched, &mut dummy.T);
            dummy.T.T = matched;
            // assert_eq!(dummy.matches(), matched * 4);
            assert_eq!(dummy.mismatches(), mismatched * 12);
        }
    }

    #[test]
    fn from_counts() {
        let mut mismatches = NucMismatches::zeros();
        mismatches.increment(
            &[Nucleotide::A, Nucleotide::Unknown, Nucleotide::C],
            &[
                NucCounts { A: 10, C: 0, G: 15, T: 0 },
                NucCounts { A: 10, C: 125, G: 15, T: 10 },
                NucCounts { A: 9, C: 0, G: 8, T: 0 },
            ],
        );
        let mut expected = NucMismatches::zeros();
        expected.A.A = 10;
        expected.A.G = 15;
        expected.C.A = 9;
        expected.C.G = 8;

        assert_eq!(mismatches, expected);
    }
}
