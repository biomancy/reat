use bio_types::strand::ReqStrand;
#[cfg(test)]
use mockall::{mock, predicate::*};
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::Record;

use crate::rada::modules::dna::Nucleotide;

pub trait Read {
    fn name(&self) -> &[u8];
    fn mapq(&self) -> u8;
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

pub trait AlignedRead: Read {
    fn cigar(&self) -> CigarStringView;
    fn pos(&self) -> i64;
    fn contig(&self) -> &str;
}

mock! {
    pub Read {}
    impl AlignedRead for Read {
        fn cigar(&self) -> CigarStringView;
        fn pos(&self) -> i64;
        fn contig(&self) -> &str;
    }

    impl Read for Read {
        fn name(&self) -> &[u8];
        fn mapq(&self) -> u8;
        fn strand(&self) -> &ReqStrand;

        fn seq(&self) -> Vec<u8>;
        fn base(&self, i: usize) -> u8;

        fn qual(&self) -> &[u8];
        fn base_qual(&self, i: usize) -> u8;

        fn is_first(&self) -> bool;
        fn len(&self) -> usize;
    }
}

impl Read for Record {
    fn name(&self) -> &[u8] {
        self.name()
    }

    fn mapq(&self) -> u8 {
        self.mapq()
    }

    fn strand(&self) -> &ReqStrand {
        if self.is_reverse() {
            &ReqStrand::Reverse
        } else {
            &ReqStrand::Forward
        }
    }

    fn seq(&self) -> Vec<u8> {
        self.seq().as_bytes()
    }

    fn qual(&self) -> &[u8] {
        self.qual()
    }

    fn is_first(&self) -> bool {
        self.is_first_in_template()
    }

    fn len(&self) -> usize {
        self.len()
    }
}

impl AlignedRead for Record {
    fn cigar(&self) -> CigarStringView {
        self.cigar()
    }

    fn pos(&self) -> i64 {
        self.pos()
    }

    fn contig(&self) -> &str {
        self.contig()
    }
}
