use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Locus, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::izip;
use itertools::multiunzip;
use itertools::Itertools;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::site::{BinnedSiteMismatches, REATBatchedSiteMismatches};
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, AggregatedNucCountsItem, NucCounter};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::SiteWorkload;

use super::base::BaseNucCounter;

pub struct IntervalNucCounts<'a> {
    bin: Interval,
    strand: Strand,
    intervals: Vec<Range<Position>>,
    seqnuc: Vec<&'a [NucCounts]>,
}

impl<'a> AggregatedNucCounts<'a> for IntervalNucCounts<'a> {
    type Meta = ();

    fn items(&'a self) -> &'a [AggregatedNucCountsItem<'a, Self::Meta>] {
        todo!()
    }

    fn consume(self) -> Vec<AggregatedNucCountsItem<'a, Self::Meta>> {
        todo!()
    }
}

impl<'a> AbstractInterval for IntervalNucCounts<'a> {
    fn contig(&self) -> &str {
        self.bin.contig()
    }

    fn range(&self) -> Range<Position> {
        self.bin.range()
    }
}

// impl<'a> AggregatedNucCounts<'a> for BinnedIntervalNucCounts<'a> {
//     fn elemrng(&self) -> &[Range<Position>] {
//         &self.intervals
//     }
//
//     fn seqnuc(&'a self) -> &'a [&'a [NucCounts]] {
//         &self.seqnuc
//     }
// }

// impl<'a> ToMismatches<'a, REATBinnedSiteMismatches> for BinnedIntervalNucCounts<'a> {
//     type MismatchesPreview = (&'a Nucleotide, &'a NucCounts);
//
//     fn mismatches(
//         self,
//         refnuc: Vec<&'a [Nucleotide]>,
//         prednuc: Vec<&'a [Nucleotide]>,
//         prefilter: impl Fn(Self::MismatchesPreview) -> bool,
//     ) -> Vec<REATBinnedSiteMismatches> {
//         debug_assert!([refnuc.len(), prednuc.len(), self.seqnuc.len(), self.intervals.len()].iter().all_equal());
//
//         let res = izip!(self.intervals.into_iter(), refnuc.into_iter(), prednuc.into_iter(), self.seqnuc.into_iter())
//             .into_iter()
//             .map(|(ranges, refn, predn, seqn)| {
//                 izip!(ranges, refn, predn, seqn).filter(|(loci, rn, pn, sn)| prefilter((pn, sn)))
//             })
//             .flatten();
//         let (loci, refnuc, prednuc, seqnuc) = multiunzip(res);
//
//         vec![REATBinnedSiteMismatches::new(self.bin.contig().into(), self.strand, loci, refnuc, prednuc, seqnuc)]
//     }
// }

pub struct IntervalNucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    ranges: Vec<Range<u64>>,
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for IntervalNucCounter<R, Filter>
where
    Self: 'a,
{
    type ColliderResult = IntervalNucCounts<'a>;
    type Workload = SiteWorkload;

    fn reset(&mut self, info: Self::Workload) {
        let (interval, ranges) = info.dissolve();
        self.base.reset(interval);
        self.ranges = ranges;
    }

    fn collide(&mut self, read: &R) {
        self.base.count(read);
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Self::ColliderResult {
        // let (contig, start) = (self.base.interval().contig(), self.base.interval().range().start);
        // let items = self
        //     .ranges
        //     .iter()
        //     .map(|x| {
        //         let indx = (x.start - start) as usize..(x.end - start) as usize;
        //         IntervalNucCounts {
        //             block: Interval::new(contig.into(), x.clone()),
        //             ncounts: &self.base.counted()[indx],
        //             strand: Strand::Unknown,
        //         }
        //     })
        //     .collect();
        // ASDAD::new(self.base.interval().clone(), items)
        todo!()
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R> for IntervalNucCounter<R, Filter> {}
