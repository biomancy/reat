use std::marker::PhantomData;
use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Interval};
use itertools::{zip, Itertools};

use crate::core::counting::buffers::{IntervalCounts, RawCounts};
use crate::core::workload::{ROIWorkload, ROI};

use super::CountsBuffer;
use super::NucCounts;

#[derive(Clone)]
pub struct ROIBuffer {
    inner: Vec<NucCounts>,
    interval: Interval,
    rois: Vec<ROI>,
    roisbuf: Vec<Interval>,
    coverage: Vec<u32>,
    index: IntervalTree<u32, usize>,
}

impl ROIBuffer {
    pub fn new(maxsize: usize) -> Self {
        let mut inner = Vec::with_capacity(maxsize);
        inner.resize(maxsize, NucCounts::zeros());
        let interval = Interval::new(String::default(), 0..maxsize as u64);
        Self { inner, interval, rois: vec![], roisbuf: vec![], coverage: vec![], index: Default::default() }
    }
}

impl CountsBuffer for ROIBuffer {
    #[inline]
    fn interval(&self) -> &Interval {
        &self.interval
    }

    fn rois(&self) -> &[Interval] {
        &self.roisbuf
    }

    #[inline]
    fn buffer(&self) -> &[NucCounts] {
        &self.inner
    }

    #[inline]
    fn buffer_mut(&mut self) -> &mut [NucCounts] {
        &mut self.inner
    }

    #[inline]
    fn add_matched(&mut self, blocks: &[Range<u32>]) {
        let index = &self.index;
        let coverage = &mut self.coverage;
        for unique in blocks.iter().map(|x| index.find(x)).flatten().map(|x| x.data()).unique() {
            coverage[*unique] += 1
        }
    }

    fn results(&self) -> Vec<IntervalCounts> {
        let mut results = Vec::with_capacity(self.rois.len());
        let interstart = self.interval.range().start;
        for (roi, cov) in zip(&self.rois, &self.coverage) {
            let (start, end) = (roi.interval.range().start - interstart, roi.interval.range().end - interstart);
            let (start, end) = (start as usize, end as usize);
            results.push(IntervalCounts {
                roi: &roi.interval,
                name: &roi.name,
                cnts: RawCounts { nuc: &self.inner[start..end], coverage: *cov },
            })
        }
        results
    }

    fn reset(&mut self, info: ROIWorkload) {
        let info = info.dissolve();
        self.interval = info.0;
        self.rois = info.1;
        self.roisbuf = self.rois.iter().map(|x| x.interval.clone()).collect();

        self.coverage.clear();
        self.coverage.resize(self.rois.len(), 0);

        let newlen = (self.interval.range().end - self.interval.range().start) as usize;
        self.inner.clear();
        self.inner.resize(newlen, NucCounts::zeros());

        // Build index
        let interstart = self.interval.range().start;
        self.index = Default::default();
        for (ind, roi) in self.rois.iter().enumerate() {
            let (start, end) = (roi.interval.range().start - interstart, roi.interval.range().end - interstart);
            let (start, end) = (start as u32, end as u32);
            self.index.insert(start..end, ind)
        }
        // self.index.index();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reset() {
        let mut dummy = ROIBuffer::new(10);
        let reset = ROIWorkload::new(
            Interval::new("".into(), 0..123),
            vec![ROI { interval: Interval::new("".into(), 10..15), name: "1".into() }],
        );

        for x in [20, 10, 5] {
            dummy.reset(reset.clone());
            // previous changes must be cleaned
            assert!(dummy.buffer().iter().all(|x| x.coverage() == 0), "{:?}", dummy.buffer());
            // new dummy changes
            dummy.buffer_mut()[0].T = 100;
        }
    }

    #[test]
    fn content() {
        let mut dummy = ROIBuffer::new(0);

        let interval = |x| Interval::new("".into(), x);
        let workload = ROIWorkload::new(
            interval(13..23),
            vec![
                ROI { interval: interval(13..16), name: "".into() }, // 0 -> 3
                ROI { interval: interval(14..18), name: "".into() }, // 1 -> 5
                ROI { interval: interval(18..22), name: "".into() }, // 5 -> 9
            ],
        );
        dummy.reset(workload.clone());

        dummy.buffer_mut()[0] = NucCounts::A(1); // 13
        dummy.buffer_mut()[1] = NucCounts::C(3); // 14
        dummy.buffer_mut()[2] = NucCounts::G(1); // 15
        dummy.buffer_mut()[3] = NucCounts::T(1); // 16
        dummy.buffer_mut()[4] = NucCounts::A(10); // 17
        dummy.buffer_mut()[5] = NucCounts::C(10); // 18
        dummy.buffer_mut()[9] = NucCounts::G(10); // 22

        dummy.add_matched(&[0..1]);
        dummy.add_matched(&[3..4, 4..5]);
        dummy.add_matched(&[5..6, 8..10]);

        dummy.add_matched(&[3..4, 5..8]);
        dummy.add_matched(&[5..6, 7..8, 9..10]);

        dummy.add_matched(&[0..2, 8..9]);
        dummy.add_matched(&[2..6]);

        let mut expected = Vec::new();
        for (ind, cov, range) in [(0usize, 3u32, 0usize..3usize), (1, 4, 1..5), (2, 5, 5..9)] {
            expected.push(IntervalCounts {
                roi: &workload.rois()[ind].interval,
                name: &workload.rois()[ind].name,
                cnts: RawCounts { nuc: &dummy.inner[range], coverage: cov },
            });
        }

        let content = dummy.results();
        assert_eq!(content, expected);
    }
}
