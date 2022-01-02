use bio_types::genome::{AbstractInterval, Interval};
use bio_types::strand::{ReqStrand, Strand};
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::interval::owned::OwnedIntervalMismatches;
use crate::core::mismatches::interval::IntervalMismatches;
use crate::core::mismatches::IntermediateMismatches;

#[derive(Constructor)]
pub struct BorrowedIntervalMismatches<'a> {
    interval: Interval,
    strand: Strand,
    refnuc: &'a [Nucleotide],
    ncounts: &'a [NucCounts],
}

impl IntervalMismatches for BorrowedIntervalMismatches<'_> {
    fn interval(&self) -> &Interval {
        &self.interval
    }

    fn strand(&self) -> &Strand {
        &self.strand
    }

    fn refnuc(&self) -> &[Nucleotide] {
        self.refnuc
    }

    fn ncounts(&self) -> &[NucCounts] {
        self.ncounts
    }

    fn split(self, indices: &[usize]) -> Vec<Self>
    where
        Self: Sized,
    {
        if indices.is_empty() {
            return vec![self];
        }

        debug_assert!(indices[0] > 0);
        let mut result = Vec::with_capacity(indices.len());

        let mut curind = 0;
        let (contig, start) = (self.interval.contig(), self.interval.range().start as u64);
        for &nextind in indices {
            result.push(Self {
                interval: Interval::new(contig.to_owned(), start + curind as u64..start + nextind as u64),
                strand: self.strand,
                refnuc: &self.refnuc[curind..nextind],
                ncounts: &self.ncounts[curind..nextind],
            });
            curind = nextind;
        }
        if *indices.last().unwrap() != self.refnuc.len() {
            let nextind = self.refnuc.len();
            result.push(Self {
                interval: Interval::new(contig.to_owned(), start + curind as u64..start + nextind as u64),
                strand: self.strand,
                refnuc: &self.refnuc[curind..],
                ncounts: &self.ncounts[curind..],
            });
        }

        debug_assert_eq!(result.len(), indices.len() + 1);
        result
    }
}

impl IntermediateMismatches for BorrowedIntervalMismatches<'_> {
    type Finalized = OwnedIntervalMismatches;

    fn finalize(self) -> Self::Finalized {
        self.into()
    }

    fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }
}
