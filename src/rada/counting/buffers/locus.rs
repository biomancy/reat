use std::cmp::max;

use derive_more::{Add, AddAssign, Constructor, Mul};

use crate::rada::dna::ReqNucleotide;

#[derive(Clone, Copy, Eq, PartialEq, Debug, Add, AddAssign, Mul, Constructor)]
#[allow(non_snake_case)]
pub struct LocusCounts {
    pub A: u32,
    pub C: u32,
    pub G: u32,
    pub T: u32,
}

impl LocusCounts {
    #[inline]
    pub fn zeros() -> LocusCounts {
        LocusCounts { A: 0, T: 0, G: 0, C: 0 }
    }

    #[inline]
    pub fn coverage(&self) -> u32 {
        self.A + self.T + self.G + self.C
    }

    #[inline]
    pub fn mismatches(&self, reference: &ReqNucleotide) -> u32 {
        match reference {
            ReqNucleotide::A => self.T + self.G + self.C,
            ReqNucleotide::C => self.A + self.T + self.G,
            ReqNucleotide::G => self.A + self.T + self.C,
            ReqNucleotide::T => self.A + self.G + self.C,
        }
    }

    #[inline]
    pub fn mostfreq(&self) -> (ReqNucleotide, &u32) {
        let maximum = max(max(self.A, self.T), max(self.G, self.C));

        if self.A == maximum {
            (ReqNucleotide::A, &self.A)
        } else if self.C == maximum {
            (ReqNucleotide::C, &self.C)
        } else if self.G == maximum {
            (ReqNucleotide::G, &self.G)
        } else {
            (ReqNucleotide::T, &self.T)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage() {
        let dummy = LocusCounts { A: 1, C: 2, G: 3, T: 0 };
        assert_eq!(dummy.coverage(), 6);
        assert_eq!(LocusCounts::zeros().coverage(), 0);
    }

    #[test]
    fn mismatches() {
        let dummy = LocusCounts { A: 1, C: 2, G: 3, T: 4 };
        assert_eq!(dummy.mismatches(&ReqNucleotide::A), 9);
        assert_eq!(dummy.mismatches(&ReqNucleotide::C), 8);
        assert_eq!(dummy.mismatches(&ReqNucleotide::G), 7);
        assert_eq!(dummy.mismatches(&ReqNucleotide::T), 6);
    }

    #[test]
    fn mostfreq_maximum() {
        let mut dummy = LocusCounts { A: 10, C: 2, G: 3, T: 5 };
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::A, &10));
        dummy.A = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::T, &5));
        dummy.T = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::G, &3));
        dummy.G = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::C, &2));
        dummy.C = 1;
    }

    #[test]
    fn mostfreq_compet_maximum() {
        let mut dummy = LocusCounts { A: 1, C: 1, G: 1, T: 1 };
        // ordered when buffers are equal
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::A, &1));
        dummy.A = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::C, &1));
        dummy.C = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::G, &1));
        dummy.G = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::T, &1));
    }

    #[test]
    fn add() {
        let mut a = LocusCounts { A: 0, C: 1, G: 2, T: 3 };
        let b = LocusCounts { A: 1, C: 2, G: 3, T: 4 };
        let result = LocusCounts { A: 1, C: 3, G: 5, T: 7 };
        assert_eq!(a + b, result);
        a += b;
        assert_eq!(a, result);
    }
}