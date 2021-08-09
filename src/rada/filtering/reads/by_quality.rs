use derive_more::Constructor;

use super::{AlignedRead, ReadsFilter};

#[derive(Constructor)]
pub struct ReadsFilterByQuality {
    mapq: u8,
    phread: u8,
}

impl<R: AlignedRead> ReadsFilter<R> for ReadsFilterByQuality {
    // 255 means that mapping quality is not available
    #[inline]
    fn is_read_ok(&self, record: &R) -> bool {
        record.mapq() != 255 && record.mapq() >= self.mapq
    }

    #[inline]
    fn is_base_ok(&self, record: &R, base: usize) -> bool {
        record.base_qual(base) >= self.phread
    }
}

#[cfg(test)]
mod tests {
    use mockall::predicate::eq;

    use crate::rada::read::MockRead;

    use super::*;

    #[test]
    fn is_read_ok() {
        let dummy = ReadsFilterByQuality::new(10, 25);

        let mut read = MockRead::new();
        for mapq in [0, 9, 255] {
            read.expect_mapq().return_const(mapq);
            assert!(!ReadsFilter::<MockRead>::is_read_ok(&dummy, &read));
            read.checkpoint();
        }
        for mapq in [10, 30, 254] {
            read.expect_mapq().return_const(mapq);
            assert!(ReadsFilter::<MockRead>::is_read_ok(&dummy, &read));
            read.checkpoint();
        }
    }

    #[test]
    fn is_base_ok() {
        let dummy = ReadsFilterByQuality::new(10, 25);

        let mut read = MockRead::new();
        for phread in [0, 24] {
            read.expect_base_qual().with(eq(1)).once().return_const(phread);
            assert!(!ReadsFilter::<MockRead>::is_base_ok(&dummy, &read, 1));
        }
        for phread in [25, 50, 255] {
            read.expect_base_qual().with(eq(1)).once().return_const(phread);
            assert!(ReadsFilter::<MockRead>::is_base_ok(&dummy, &read, 1));
        }
    }
}
