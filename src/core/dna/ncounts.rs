use std::cmp::Ordering;
use std::ops::{Index, IndexMut};

use derive_more::{Add, AddAssign, Mul};

use crate::core::dna::{Nucleotide, ReqNucleotide};
use funty::Numeric;

pub type NucCounts = InnerNucCounts<u32>;
pub type FracNucCounts = InnerNucCounts<f32>;

#[derive(Clone, Copy, Eq, PartialEq, Debug, Add, AddAssign, Mul, Default)]
#[allow(non_snake_case)]
pub struct InnerNucCounts<T: Numeric + Sized> {
    pub A: T,
    pub C: T,
    pub G: T,
    pub T: T,
}

impl<Num: Numeric> InnerNucCounts<Num> {
    #[allow(non_snake_case)]
    pub fn new(A: Num, C: Num, G: Num, T: Num) -> Self {
        Self { A, C, G, T }
    }

    // pub fn increment(&mut self, sequence: &[Nucleotide]) {
    //     for nuc in sequence {
    //         match nuc {
    //             Nucleotide::A => self.A += 1.into(),
    //             Nucleotide::C => self.C += 1.into(),
    //             Nucleotide::G => self.G += 1.into(),
    //             Nucleotide::T => self.T += 1.into(),
    //             Nucleotide::Unknown => {}
    //         }
    //     }
    // }

    #[inline]
    pub fn zeros() -> Self {
        Self { A: Num::default(), T: Num::default(), G: Num::default(), C: Num::default() }
    }

    #[inline]
    #[allow(non_snake_case)]
    pub fn A(A: Num) -> Self {
        Self { A, T: Num::default(), G: Num::default(), C: Num::default() }
    }

    #[inline]
    #[allow(non_snake_case)]
    pub fn T(T: Num) -> Self {
        Self { A: Num::default(), T, G: Num::default(), C: Num::default() }
    }

    #[inline]
    #[allow(non_snake_case)]
    pub fn G(G: Num) -> Self {
        Self { A: Num::default(), T: Num::default(), G, C: Num::default() }
    }

    #[inline]
    #[allow(non_snake_case)]
    pub fn C(C: Num) -> Self {
        Self { A: Num::default(), T: Num::default(), G: Num::default(), C }
    }

    #[inline]
    pub fn coverage(&self) -> Num {
        self.A + self.T + self.G + self.C
    }

    #[inline]
    pub fn mismatches(&self, reference: Nucleotide) -> Num {
        match reference {
            Nucleotide::A => self.C + self.G + self.T,
            Nucleotide::C => self.A + self.G + self.T,
            Nucleotide::G => self.A + self.C + self.T,
            Nucleotide::T => self.A + self.C + self.G,
            Nucleotide::Unknown => self.A + self.C + self.G + self.T,
        }
    }

    #[inline]
    pub fn mostfreq(&self) -> (ReqNucleotide, &Num) {
        let ac = match self.A.partial_cmp(&self.C) {
            Some(Ordering::Equal | Ordering::Greater) => Some(ReqNucleotide::A),
            Some(Ordering::Less) => Some(ReqNucleotide::C),
            _ => None,
        };
        let gt = match self.G.partial_cmp(&self.T) {
            Some(Ordering::Equal | Ordering::Greater) => Some(ReqNucleotide::G),
            Some(Ordering::Less) => Some(ReqNucleotide::T),
            _ => None,
        };
        if ac.is_none() && gt.is_none() {
            return (ReqNucleotide::A, &self.A);
        }
        let (ac, gt) = (ac.unwrap(), gt.unwrap());
        match self[ac].partial_cmp(&self[gt]) {
            None => (ReqNucleotide::A, &self.A),
            Some(Ordering::Equal | Ordering::Greater) => (ac, &self[ac]),
            Some(Ordering::Less) => (gt, &self[gt]),
        }
    }

    #[inline]
    pub fn complementary(&self) -> Self {
        Self { A: self.T, C: self.G, G: self.C, T: self.A }
    }
}

impl<T: Numeric> Index<ReqNucleotide> for InnerNucCounts<T> {
    type Output = T;

    fn index(&self, index: ReqNucleotide) -> &Self::Output {
        match index {
            ReqNucleotide::A => &self.A,
            ReqNucleotide::C => &self.C,
            ReqNucleotide::G => &self.G,
            ReqNucleotide::T => &self.T,
        }
    }
}

impl<T: Numeric> IndexMut<ReqNucleotide> for InnerNucCounts<T> {
    fn index_mut(&mut self, index: ReqNucleotide) -> &mut Self::Output {
        match index {
            ReqNucleotide::A => &mut self.A,
            ReqNucleotide::C => &mut self.C,
            ReqNucleotide::G => &mut self.G,
            ReqNucleotide::T => &mut self.T,
        }
    }
}

impl From<&'_ NucCounts> for FracNucCounts {
    fn from(nc: &'_ NucCounts) -> Self {
        Self { A: nc.A as f32, C: nc.C as f32, G: nc.G as f32, T: nc.T as f32 }
    }
}

impl From<NucCounts> for FracNucCounts {
    fn from(nc: NucCounts) -> Self {
        Self { A: nc.A as f32, C: nc.C as f32, G: nc.G as f32, T: nc.T as f32 }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::dna::{Nucleotide, ReqNucleotide};

    use super::*;

    // #[test]
    // fn from_sequence() {
    //     let mut counts = InnerNucCounts { A: 1, C: 2, G: 3, T: 5 };
    //
    //     let sequence = [
    //         [Nucleotide::A].repeat(7),
    //         [Nucleotide::C].repeat(10),
    //         [Nucleotide::Unknown].repeat(5),
    //         [Nucleotide::G].repeat(3),
    //         [Nucleotide::T].repeat(2),
    //         [Nucleotide::Unknown].repeat(5),
    //     ]
    //     .iter()
    //     .flatten()
    //     .copied()
    //     .collect_vec();
    //     counts.increment(&sequence);
    //
    //     let expected = InnerNucCounts { A: 8, C: 12, G: 6, T: 7 };
    //     assert_eq!(expected, counts);
    // }

    #[test]
    fn coverage() {
        let dummy = InnerNucCounts { A: 1, C: 2, G: 3, T: 0 };
        assert_eq!(dummy.coverage(), 6);
        assert_eq!(InnerNucCounts::<u32>::zeros().coverage(), 0);
    }

    #[test]
    fn mismatches() {
        let dummy = InnerNucCounts { A: 1, C: 2, G: 3, T: 4 };
        assert_eq!(dummy.mismatches(Nucleotide::A), 9);
        assert_eq!(dummy.mismatches(Nucleotide::C), 8);
        assert_eq!(dummy.mismatches(Nucleotide::G), 7);
        assert_eq!(dummy.mismatches(Nucleotide::T), 6);
        assert_eq!(dummy.mismatches(Nucleotide::Unknown), 10);
    }

    #[test]
    fn mostfreq_maximum() {
        let mut dummy = InnerNucCounts { A: 10, C: 2, G: 3, T: 5 };
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
        let mut dummy = InnerNucCounts { A: 1, C: 1, G: 1, T: 1 };
        // ordered when ncounters are equal
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
        let mut a = InnerNucCounts { A: 0, C: 1, G: 2, T: 3 };
        let b = InnerNucCounts { A: 1, C: 2, G: 3, T: 4 };
        let result = InnerNucCounts { A: 1, C: 3, G: 5, T: 7 };
        assert_eq!(a + b, result);
        a += b;
        assert_eq!(a, result);
    }
}
