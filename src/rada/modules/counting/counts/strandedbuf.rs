use bio_types::strand::ReqStrand;
use rust_htslib::bam::record::Record;

use crate::rada::modules::stranding::deduct::StrandDeductor;

use super::{CountsBuffer, CountsBufferContent, LocusCounts};
use crate::rada::modules::read::{AlignedRead, MockRead};
use std::marker::PhantomData;

pub struct StrandedCountsBuffer<R: AlignedRead, Deductor: StrandDeductor<R>> {
    forward: Vec<LocusCounts>,
    reverse: Vec<LocusCounts>,
    strand_deductor: Deductor,
    phantom: PhantomData<R>,
}

impl<R: AlignedRead, Deductor: StrandDeductor<R>> StrandedCountsBuffer<R, Deductor> {
    fn new(maxsize: u32, strand_deductor: Deductor) -> Self {
        let (mut plstrand, mut mnstrand) = (Vec::new(), Vec::new());
        plstrand.reserve(maxsize as usize);
        mnstrand.reserve(maxsize as usize);
        StrandedCountsBuffer {
            forward: plstrand,
            reverse: mnstrand,
            strand_deductor,
            phantom: Default::default(),
        }
    }
}

impl<R: AlignedRead, Deductor: StrandDeductor<R>> CountsBuffer<R>
    for StrandedCountsBuffer<R, Deductor>
{
    fn reset(&mut self, newlen: u32) -> &mut Self {
        let newlen = newlen as usize;
        if self.forward.len() != newlen {
            self.forward.resize(newlen, LocusCounts::zeros());
            self.reverse.resize(newlen, LocusCounts::zeros());
        }
        self.forward.fill(LocusCounts::zeros());
        self.reverse.fill(LocusCounts::zeros());
        self
    }

    #[inline]
    fn buffer_for(&mut self, record: &mut R) -> &mut [LocusCounts] {
        match self.strand_deductor.deduce(record) {
            ReqStrand::Forward => &mut self.forward,
            ReqStrand::Reverse => &mut self.reverse,
        }
    }

    fn content(&self) -> CountsBufferContent {
        CountsBufferContent {
            forward: Some(&self.forward),
            reverse: Some(&self.reverse),
            unstranded: None,
        }
    }

    fn len(&self) -> u32 {
        self.forward.len() as u32
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use crate::rada::modules::stranding::deduct::MockStrandDeductor;

    use super::*;

    fn dummy(
        maxsize: u32,
        deductor: MockStrandDeductor<MockRead>,
    ) -> StrandedCountsBuffer<MockRead, MockStrandDeductor<MockRead>> {
        StrandedCountsBuffer::<MockRead, MockStrandDeductor<MockRead>>::new(
            maxsize,
            MockStrandDeductor::new(),
        )
    }

    #[test]
    fn reset() {
        let mut dummy = dummy(10, Default::default());
        assert_eq!(dummy.len(), 0);
        for x in [20, 10, 5] {
            dummy.reset(x);
            assert_eq!(dummy.len(), x);
            // previous changes must be cleaned
            assert!(dummy.forward.iter().all(|x| x.coverage() == 0));
            assert!(dummy.reverse.iter().all(|x| x.coverage() == 0));
            // new dummy changes
            dummy.forward[0].T = 100;
            dummy.reverse.last_mut().unwrap().T = 100;
        }
    }

    #[test]
    fn buffer_for() {
        let mut mock = MockStrandDeductor::new();
        mock.expect_deduce().once().return_const(ReqStrand::Forward);
        mock.expect_deduce().once().return_const(ReqStrand::Reverse);

        let mut dummy = StrandedCountsBuffer::new(10, mock);

        let mut record = Record::new();
        assert!(ptr::eq(
            dummy.buffer_for(&mut record),
            dummy.forward.as_slice()
        ));
        assert!(ptr::eq(
            dummy.buffer_for(&mut record),
            dummy.reverse.as_slice()
        ));
    }

    #[test]
    fn content() {
        let mut dummy = dummy(10, Default::default());
        dummy.reset(10);
        dummy.forward[0].A = 10;
        dummy.reverse.last_mut().unwrap().C = 1;

        let content = dummy.content();
        assert!(content.unstranded.is_none());
        assert!(
            content.forward.is_some()
                && content.forward.unwrap().len() == 10
                && content.forward.unwrap()[0].A == 10
        );
        assert!(
            content.reverse.is_some()
                && content.reverse.unwrap().len() == 10
                && content.reverse.unwrap().last().unwrap().C == 1
        );
    }
}
