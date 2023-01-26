use std::ops::{Index, IndexMut};

use derive_more::{Add, AddAssign};

pub use inner::ROINucCounts;

use crate::core::dna::FracNucCounts;
use crate::core::dna::{Nucleotide, ReqNucleotide};

// Simple struct-level attribute is not working for some reason
mod inner {
    #![allow(non_snake_case)]

    use super::*;

    // All counts are doubled -> 2 for homozygous site / 1 for heterozygous
    #[derive(Copy, Clone, PartialEq, Debug, Add, AddAssign)]
    pub struct ROINucCounts {
        pub A: FracNucCounts,
        pub C: FracNucCounts,
        pub G: FracNucCounts,
        pub T: FracNucCounts,
    }
}

impl ROINucCounts {
    #[inline]
    pub fn zeros() -> Self {
        ROINucCounts {
            A: FracNucCounts::zeros(),
            C: FracNucCounts::zeros(),
            G: FracNucCounts::zeros(),
            T: FracNucCounts::zeros(),
        }
    }

    #[inline]
    pub fn coverage(&self) -> f32 {
        self.A.coverage() + self.C.coverage() + self.G.coverage() + self.T.coverage()
    }

    #[inline]
    pub fn mismatches(&self) -> f32 {
        self.A.mismatches(Nucleotide::A)
            + self.C.mismatches(Nucleotide::C)
            + self.G.mismatches(Nucleotide::G)
            + self.T.mismatches(Nucleotide::T)
    }

    #[inline]
    pub fn complementary(&self) -> Self {
        ROINucCounts {
            A: self.T.complementary(),
            C: self.G.complementary(),
            G: self.C.complementary(),
            T: self.A.complementary(),
        }
    }
}

impl Default for ROINucCounts {
    fn default() -> Self {
        ROINucCounts::zeros()
    }
}

impl Index<ReqNucleotide> for ROINucCounts {
    type Output = FracNucCounts;

    fn index(&self, index: ReqNucleotide) -> &Self::Output {
        match index {
            ReqNucleotide::A => &self.A,
            ReqNucleotide::C => &self.C,
            ReqNucleotide::G => &self.G,
            ReqNucleotide::T => &self.T,
        }
    }
}

impl IndexMut<ReqNucleotide> for ROINucCounts {
    fn index_mut(&mut self, index: ReqNucleotide) -> &mut Self::Output {
        match index {
            ReqNucleotide::A => &mut self.A,
            ReqNucleotide::C => &mut self.C,
            ReqNucleotide::G => &mut self.G,
            ReqNucleotide::T => &mut self.T,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fillall(value: f32, counts: &mut FracNucCounts) {
        counts.A = value;
        counts.C = value;
        counts.G = value;
        counts.T = value;
    }

    #[test]
    fn coverage() {
        let mut dummy: ROINucCounts = Default::default();
        for num in [0_f32, 25_f32] {
            fillall(num, &mut dummy.A);
            fillall(num, &mut dummy.C);
            fillall(num, &mut dummy.G);
            fillall(num, &mut dummy.T);
            assert_eq!(dummy.coverage(), num * 4_f32 * 4_f32);
        }
    }

    #[test]
    fn match_mismatch() {
        let mut dummy: ROINucCounts = Default::default();
        for (mismatched, matched) in [(0_f32, 25_f32), (12_f32, 0_f32), (1_f32, 2_f32)] {
            fillall(mismatched, &mut dummy.A);
            dummy.A.A = matched;
            fillall(mismatched, &mut dummy.C);
            dummy.C.C = matched;
            fillall(mismatched, &mut dummy.G);
            dummy.G.G = matched;
            fillall(mismatched, &mut dummy.T);
            dummy.T.T = matched;
            // assert_eq!(dummy.matches(), matched * 4);
            assert_eq!(dummy.mismatches(), mismatched * 12_f32);
        }
    }

    // #[test]
    // fn from_counts() {
    //     let mut mismatches = NucMismatches::zeros();
    //     mismatches.increment(
    //         &[Nucleotide::A, Nucleotide::Unknown, Nucleotide::C],
    //         &[
    //             NucCounts { A: 10, C: 0, G: 15, T: 0 },
    //             NucCounts { A: 10, C: 125, G: 15, T: 10 },
    //             NucCounts { A: 9, C: 0, G: 8, T: 0 },
    //         ],
    //     );
    //     let mut expected = NucMismatches::zeros();
    //     expected.A.A = 10;
    //     expected.A.G = 15;
    //     expected.C.A = 9;
    //     expected.C.G = 8;
    //
    //     assert_eq!(mismatches, expected);
    // }
}
