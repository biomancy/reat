use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Interval, Locus, Position};
use bio_types::strand::{ReqStrand, Strand};
use itertools::{zip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{BinnedROIMismatches, REATBinnedROIMismatches};
use crate::core::read::AlignedRead;
use crate::core::rpileup::ncounters::filters::ReadsFilter;
use crate::core::rpileup::ncounters::{AggregatedNucCounts, AggregatedNucCountsItem, NucCounter};
use crate::core::rpileup::ReadsCollider;
use crate::core::workload::{ROIWorkload, ROI};

use super::base::BaseNucCounter;

pub struct BinnedROINucCounts<'a> {
    bin: Interval,
    strand: Strand,
    coverage: Vec<u32>,
    // Flattened ROI structure
    rois: Vec<Range<Position>>,
    pieces: Vec<Vec<Range<Position>>>,
    premasked: Vec<Range<Position>>,
    names: Vec<String>,
    // Counts in each roi
    seqnuc: Vec<&'a [NucCounts]>,
}

impl<'a> AbstractInterval for BinnedROINucCounts<'a> {
    fn contig(&self) -> &str {
        self.bin.contig()
    }

    fn range(&self) -> Range<Position> {
        self.bin.range()
    }
}

impl<'a> AggregatedNucCounts<'a> for BinnedROINucCounts<'a> {
    type Meta = (u32, ROI);

    fn items(&'a self) -> &'a [AggregatedNucCountsItem<'a, Self::Meta>] {
        todo!()
    }

    fn consume(self) -> Vec<AggregatedNucCountsItem<'a, Self::Meta>> {
        todo!()
    }
}

// impl<'a> AggregatedNucCounts<'a> for BinnedROINucCounts<'a> {
//     fn elemrng(&self) -> &[Range<Position>] {
//         &self.rois
//     }
//
//     fn seqnuc(&'a self) -> &'a [&'a [NucCounts]] {
//         &self.seqnuc
//     }
// }

// impl<'a, 'b> ToMismatches<'a, REATBinnedROIMismatches> for BinnedROINucCounts<'a> {
//     type MismatchesPreview = ();
//
//     fn mismatches(
//         self,
//         refnuc: Vec<&'a [Nucleotide]>,
//         prednuc: Vec<&'a [Nucleotide]>,
//         prefilter: impl Fn(Self::MismatchesPreview) -> bool,
//     ) -> Vec<REATBinnedROIMismatches> {
//         todo!()
//         // debug_assert!(
//         //     (roi.range().end - roi.range().start) as usize == reference.len() && reference.len() == sequenced.len()
//         // );
//         // // Calculate mismatches only over the retained subintervals
//         // let (mut cnts, mut mismatches) = (NucCounts::zeros(), MismatchesSummary::zeros());
//         // let start = roi.range().start;
//         // let mut nucin = 0;
//         // for piece in roi.include() {
//         //     nucin += piece.end - piece.start;
//         //
//         //     let idx = (piece.start - start) as usize..(piece.end - start) as usize;
//         //     cnts.increment(&reference[idx.clone()]);
//         //     mismatches.increment(&reference[idx.clone()], &sequenced[idx])
//         // }
//         //
//         // let total = roi.original().range().end - roi.original().range().start;
//         // Self { roi, strand, masked: (total - nucin) as u32, coverage, sequence: cnts, mismatches }
//     }
// }

pub struct ROINucCounter<R: AlignedRead, Filter: ReadsFilter<R>> {
    base: BaseNucCounter<R, Filter>,
    rois: Vec<ROI>,
    coverage: Vec<u32>,
    index: IntervalTree<u32, usize>,
}

impl<R: AlignedRead, Filter: ReadsFilter<R>> ROINucCounter<R, Filter> {
    pub fn new(maxbuf: usize, filter: Filter, trim5: u16, trim3: u16) -> Self {
        Self {
            base: BaseNucCounter::new(maxbuf, filter, trim5, trim3),
            rois: vec![],
            coverage: vec![],
            index: Default::default(),
        }
    }
}

impl<'a, R: AlignedRead, Filter: ReadsFilter<R>> ReadsCollider<'a, R> for ROINucCounter<R, Filter> {
    type ColliderResult = BinnedROINucCounts<'a>;
    type Workload = ROIWorkload;

    fn reset(&mut self, info: Self::Workload) {
        let (interval, rois) = info.dissolve();
        self.base.reset(interval);
        self.rois = rois;

        // Coverage for each roi
        self.coverage.clear();
        self.coverage.resize(self.rois.len(), 0);

        // Index to accurately count ROIs coverage
        self.index = Default::default();
        let binstart = self.base.interval().range().start;
        for (ind, roi) in self.rois.iter().enumerate() {
            let (start, end) = (roi.range().start - binstart, roi.range().end - binstart);
            self.index.insert(start as u32..end as u32, ind)
        }
    }

    fn collide(&mut self, read: &R) {
        let hits = self.base.count(read).iter().map(|x| self.index.find(x)).flatten().map(|x| x.data()).unique();
        for ind in hits {
            self.coverage[*ind] += 1;
        }
    }

    fn finalize(&mut self) {}

    fn result(&'a self) -> Self::ColliderResult {
        // let ncounts = self.base.counted();
        // let interval = self.base.interval();
        // let binstart = interval.range().start as usize;
        //
        // let mut results = Vec::with_capacity(self.rois.len());
        // for (coverage, roi) in zip(&self.coverage, &self.rois) {
        //     debug_assert_eq!(roi.original().contig(), interval.contig());
        //     let (start, end) = (roi.range().start as usize, roi.range().end as usize);
        //     let roicnts = &ncounts[start - binstart..end - binstart];
        //     results.push(ROINucCounts { coverage: *coverage, roi, strand: Strand::Unknown, ncounts: roicnts });
        // }
        // ASDAD::new(interval.clone(), results)
        todo!()
    }
}

impl<'a, R: 'a + AlignedRead, Filter: 'a + ReadsFilter<R>> NucCounter<'a, R> for ROINucCounter<R, Filter> {}
