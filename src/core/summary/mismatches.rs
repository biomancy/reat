pub use inner::*;

use crate::core::dna::Nucleotide;

// See NucCounts note
mod inner {
    #![allow(non_snake_case)]

    use derive_more::{Add, AddAssign, Constructor};

    pub use crate::core::counting::NucCounts;

    // These are counts with respect to the FORWARD strand
    #[derive(Constructor, Copy, Clone, Eq, PartialEq, Debug, Add, AddAssign)]
    pub struct MismatchesSummary {
        pub A: NucCounts,
        pub C: NucCounts,
        pub G: NucCounts,
        pub T: NucCounts,
    }
}

impl MismatchesSummary {
    #[inline]
    pub fn zeros() -> Self {
        MismatchesSummary { A: NucCounts::zeros(), C: NucCounts::zeros(), G: NucCounts::zeros(), T: NucCounts::zeros() }
    }

    #[inline]
    pub fn coverage(&self) -> u32 {
        self.A.coverage() + self.C.coverage() + self.G.coverage() + self.T.coverage()
    }

    #[inline]
    pub fn mismatches(&self) -> u32 {
        self.A.mismatches(&Nucleotide::A)
            + self.C.mismatches(&Nucleotide::C)
            + self.G.mismatches(&Nucleotide::G)
            + self.T.mismatches(&Nucleotide::T)
    }

    #[inline]
    pub fn complementary(&self) -> Self {
        MismatchesSummary {
            A: self.T.complementary(),
            C: self.G.complementary(),
            G: self.C.complementary(),
            T: self.A.complementary(),
        }
    }
}

impl Default for MismatchesSummary {
    fn default() -> Self {
        MismatchesSummary::zeros()
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
        let mut dummy: MismatchesSummary = Default::default();
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
        let mut dummy: MismatchesSummary = Default::default();
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
}
