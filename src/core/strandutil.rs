use std::ops::{Index, IndexMut};

use bio_types::strand::Strand;

#[derive(Default)]
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
