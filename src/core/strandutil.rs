use std::fmt::{Display, Formatter};
use std::iter::Sum;
use std::ops::{Index, IndexMut};

use bio_types::strand::Strand;
use derive_more::{Add, AddAssign};

#[derive(Default, Copy, Clone, Add, AddAssign)]
pub struct Stranded<T> {
    pub forward: T,
    pub reverse: T,
    pub unknown: T,
}

impl<T> Index<Strand> for Stranded<T> {
    type Output = T;

    fn index(&self, index: Strand) -> &Self::Output {
        match index {
            Strand::Forward => &self.forward,
            Strand::Reverse => &self.reverse,
            Strand::Unknown => &self.unknown,
        }
    }
}

impl<T> IndexMut<Strand> for Stranded<T> {
    fn index_mut(&mut self, index: Strand) -> &mut Self::Output {
        match index {
            Strand::Forward => &mut self.forward,
            Strand::Reverse => &mut self.reverse,
            Strand::Unknown => &mut self.unknown,
        }
    }
}

impl<T: Display> Display for Stranded<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}(+), {}(-), {}(.))", self.forward, self.reverse, self.unknown)
    }
}

impl<T> Stranded<T> {
    #[inline(always)]
    pub fn with_fn(func: impl Fn(Strand) -> T) -> Self {
        Self { forward: func(Strand::Forward), reverse: func(Strand::Reverse), unknown: func(Strand::Unknown) }
    }

    #[inline(always)]
    pub fn apply_mut(&mut self, func: impl Fn(&mut T, Strand)) {
        func(&mut self.forward, Strand::Forward);
        func(&mut self.reverse, Strand::Reverse);
        func(&mut self.unknown, Strand::Unknown);
    }

    #[inline(always)]
    pub fn apply(&self, func: impl Fn(&T, Strand)) {
        func(&self.forward, Strand::Forward);
        func(&self.reverse, Strand::Reverse);
        func(&self.unknown, Strand::Unknown);
    }

    #[inline(always)]
    pub fn into<I>(self, func: impl Fn(T, Strand) -> I) -> Stranded<I> {
        Stranded {
            forward: func(self.forward, Strand::Forward),
            reverse: func(self.reverse, Strand::Reverse),
            unknown: func(self.unknown, Strand::Unknown),
        }
    }
}

impl<T: Default> Stranded<T> {
    pub fn forward(data: T) -> Self {
        Self { forward: data, reverse: Default::default(), unknown: Default::default() }
    }

    pub fn reverse(data: T) -> Self {
        Self { forward: Default::default(), reverse: data, unknown: Default::default() }
    }

    pub fn unknown(data: T) -> Self {
        Self { forward: Default::default(), reverse: Default::default(), unknown: data }
    }
}
