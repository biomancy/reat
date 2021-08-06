use std::cmp::max;

use derive_more::{Add, AddAssign};

use crate::rada::modules::dna::ReqNucleotide;

#[derive(Clone, Copy, Add, AddAssign)]
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
        LocusCounts {
            A: 0,
            T: 0,
            G: 0,
            C: 0,
        }
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
    fn test_coverage() {
        let mut dummy = LocusCounts {
            A: 1,
            C: 2,
            G: 3,
            T: 0,
        };
        assert_eq!(dummy.coverage(), 6);
        assert_eq!(LocusCounts::zeros().coverage(), 0);
    }

    #[test]
    fn test_mismatches() {
        let dummy = LocusCounts {
            A: 1,
            C: 2,
            G: 3,
            T: 4,
        };
        assert_eq!(dummy.mismatches(&ReqNucleotide::A), 9);
        assert_eq!(dummy.mismatches(&ReqNucleotide::C), 8);
        assert_eq!(dummy.mismatches(&ReqNucleotide::G), 7);
        assert_eq!(dummy.mismatches(&ReqNucleotide::T), 6);
    }

    #[test]
    fn test_mostfreq() {
        let mut dummy = LocusCounts {
            A: 10,
            C: 2,
            G: 3,
            T: 5,
        };
        // naive maximum
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::A, &10));
        dummy.A = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::T, &5));
        dummy.T = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::G, &3));
        dummy.G = 1;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::C, &2));
        dummy.C = 1;
        // ordered when counts are equal
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::A, &1));
        dummy.A = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::C, &1));
        dummy.C = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::G, &1));
        dummy.G = 0;
        assert_eq!(dummy.mostfreq(), (ReqNucleotide::T, &1));
    }
}
