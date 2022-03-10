use std::ops::{Index, IndexMut};

use bio_types::strand::Strand;

#[derive(Default, Copy, Clone)]
pub struct StrandedData<T> {
    pub forward: T,
    pub reverse: T,
    pub unknown: T,
}

impl<T> Index<Strand> for StrandedData<T> {
    type Output = T;

    fn index(&self, index: Strand) -> &Self::Output {
        match index {
            Strand::Forward => &self.forward,
            Strand::Reverse => &self.reverse,
            Strand::Unknown => &self.unknown,
        }
    }
}

impl<T> IndexMut<Strand> for StrandedData<T> {
    fn index_mut(&mut self, index: Strand) -> &mut Self::Output {
        match index {
            Strand::Forward => &mut self.forward,
            Strand::Reverse => &mut self.reverse,
            Strand::Unknown => &mut self.unknown,
        }
    }
}

impl<T> StrandedData<T> {
    #[inline(always)]
    pub fn apply_mut(&mut self, func: impl Fn(T, Strand) -> T)
    where
        T: Default,
    {
        self.forward = func(std::mem::take(&mut self.forward), Strand::Forward);
        self.reverse = func(std::mem::take(&mut self.reverse), Strand::Reverse);
        self.unknown = func(std::mem::take(&mut self.unknown), Strand::Unknown);
    }

    #[inline(always)]
    pub fn apply(&self, func: impl Fn(&T, Strand)) {
        func(&self.forward, Strand::Forward);
        func(&self.reverse, Strand::Reverse);
        func(&self.unknown, Strand::Unknown);
    }
}
