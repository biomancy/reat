use bio_types::genome;
use bio_types::strand::ReqStrand;
#[cfg(test)]
use mockall::{mock, predicate::*};
use rust_htslib::bam::Record;
use rust_htslib::bam::record::CigarStringView;

#[allow(clippy::len_without_is_empty)]
pub trait SequencedRead {
    fn name(&self) -> &[u8];
    fn strand(&self) -> &ReqStrand;

    fn seq(&self) -> Vec<u8>;
    fn base(&self, i: usize) -> u8 {
        self.seq()[i]
    }

    fn qual(&self) -> &[u8];
    fn base_qual(&self, i: usize) -> u8 {
        self.qual()[i]
    }

    fn is_first(&self) -> bool;

    fn len(&self) -> usize;
}

pub trait AlignedRead: SequencedRead {
    fn cigar(&self) -> CigarStringView;
    fn mapq(&self) -> u8;
    fn pos(&self) -> i64;
    fn contig(&self) -> &str;
    fn flags(&self) -> u16;
}

#[cfg(test)]
mock! {
    pub Read {}
    impl AlignedRead for Read {
        fn cigar(&self) -> CigarStringView;
        fn mapq(&self) -> u8;
        fn pos(&self) -> i64;
        fn contig(&self) -> &str;
        fn flags(&self) -> u16;
    }

    impl SequencedRead for Read {
        fn name(&self) -> &[u8];
        fn strand(&self) -> &ReqStrand;

        fn seq(&self) -> Vec<u8>;
        fn base(&self, i: usize) -> u8;

        fn qual(&self) -> &[u8];
        fn base_qual(&self, i: usize) -> u8;

        fn is_first(&self) -> bool;
        fn len(&self) -> usize;
    }
}

impl SequencedRead for Record {
    #[inline]
    fn name(&self) -> &[u8] {
        self.qname()
    }

    #[inline]
    fn strand(&self) -> &ReqStrand {
        if self.is_reverse() {
            &ReqStrand::Reverse
        } else {
            &ReqStrand::Forward
        }
    }

    #[inline]
    fn seq(&self) -> Vec<u8> {
        self.seq().as_bytes()
    }

    #[inline]
    fn qual(&self) -> &[u8] {
        self.qual()
    }

    #[inline]
    fn is_first(&self) -> bool {
        self.is_first_in_template()
    }

    #[inline]
    fn len(&self) -> usize {
        self.seq_len()
    }
}

impl AlignedRead for Record {
    #[inline]
    fn cigar(&self) -> CigarStringView {
        self.cigar()
    }

    #[inline]
    fn mapq(&self) -> u8 {
        self.mapq()
    }
    #[inline]
    fn pos(&self) -> i64 {
        self.pos()
    }

    #[inline]
    fn contig(&self) -> &str {
        genome::AbstractInterval::contig(self)
    }

    #[inline]
    fn flags(&self) -> u16 {
        self.flags()
    }
}
