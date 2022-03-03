use derive_getters::Getters;
use derive_more::Constructor;

use super::{AlignedRead, ReadsFilter};

#[derive(Constructor, Getters, Copy, Clone)]
pub struct ByFlags {
    include: u16,
    exclude: u16,
}

impl<R: AlignedRead> ReadsFilter<R> for ByFlags {
    #[inline]
    fn is_read_ok(&self, record: &R) -> bool {
        ((record.flags() & self.include) == self.include) && ((record.flags() & self.exclude) == 0)
    }

    #[inline]
    fn is_base_ok(&self, _: &R, _: usize) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::core::read::MockRead;

    use super::*;

    #[test]
    fn is_read_ok() {
        let dummy = ByFlags::new(3, 3844);

        let mut read = MockRead::new();

        for (flag, result) in [(83u16, true), (91u16, true), (1107u16, false), (1, false), (4095, false)] {
            read.expect_flags().return_const(flag);
            assert_eq!(dummy.is_read_ok(&read), result);
            read.checkpoint()
        }
    }
}
