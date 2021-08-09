use derive_more::Constructor;

use crate::rada::counting::LocusCounts;
use crate::rada::dna::ReqNucleotide;

#[derive(Constructor)]
#[allow(non_snake_case)]
pub struct MismatchesSummary {
    pub A: LocusCounts,
    pub C: LocusCounts,
    pub G: LocusCounts,
    pub T: LocusCounts,
}

impl MismatchesSummary {
    pub fn zeros() -> Self {
        MismatchesSummary {
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

impl Default for MismatchesSummary {
    fn default() -> Self {
        MismatchesSummary::zeros()
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use bio_types::strand::Strand;
    use mockall::predicate::eq;

    use super::*;

    fn fillall(value: u32, counts: &mut LocusCounts) {
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
            assert_eq!(dummy.matches(), matched * 4);
            assert_eq!(dummy.mismatches(), mismatched * 12);
        }
    }
}
