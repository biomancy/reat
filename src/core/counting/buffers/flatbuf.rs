
use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval};

use crate::core::counting::buffers::{IntervalCounts, RawCounts};

use super::CountsBuffer;
use super::NucCounts;
use crate::core::workload::ROIWorkload;

#[derive(Clone)]
pub struct FlatBuffer {
    inner: Vec<NucCounts>,
    interval: Interval,
    rois: Vec<Interval>,
    coverage: u32,
}

impl FlatBuffer {
    pub fn new(maxsize: usize) -> Self {
        let mut inner = Vec::with_capacity(maxsize);
        inner.resize(maxsize, NucCounts::zeros());
        let interval = Interval::new(String::default(), 0..maxsize as u64);
        let rois: Vec<Interval> = vec![interval.clone()];
        Self { inner, interval, rois, coverage: 0 }
    }
}

impl CountsBuffer for FlatBuffer {
    #[inline]
    fn interval(&self) -> &Interval {
        &self.interval
    }

    fn rois(&self) -> &[Interval] {
        &self.rois
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
    fn add_matched(&mut self, _: &[Range<u32>]) {
        self.coverage += 1;
    }

    #[inline]
    fn results(&self) -> Vec<IntervalCounts> {
        vec![IntervalCounts {
            roi: &self.interval,
            name: "",
            cnts: RawCounts { nuc: &self.inner, coverage: self.coverage },
        }]
    }

    fn reset(&mut self, workload: ROIWorkload) {
        let (interval, _) = workload.dissolve();
        let newlen = interval.range().end - interval.range().start;
        debug_assert!(newlen > 0);
        self.inner.clear();
        self.inner.resize(newlen as usize, NucCounts::zeros());
        self.rois = vec![interval.clone()];
        self.interval = interval;
        self.coverage = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reset() {
        let mut dummy = FlatBuffer::new(10);
        for x in [20, 10, 5] {
            dummy.reset(ROIWorkload::new(Interval::new("".into(), 0..x), vec![]));
            // previous changes must be cleaned
            assert!(dummy.buffer().iter().all(|x| x.coverage() == 0), "{:?}", dummy.buffer());
            // new dummy changes
            dummy.buffer_mut()[0].T = 100;
        }
    }

    #[test]
    fn content() {
        let mut dummy = FlatBuffer::new(10);
        let interval = Interval::new("".into(), 0..3);
        dummy.reset(ROIWorkload::new(interval.clone(), vec![]));
        dummy.buffer_mut()[0].A = 10;
        dummy.buffer_mut()[2].C = 3;
        dummy.add_matched(&[0..1, 2..3]);

        let counts = [NucCounts { A: 10, C: 0, G: 0, T: 0 }, NucCounts::zeros(), NucCounts { A: 0, C: 3, G: 0, T: 0 }];
        let expected = vec![IntervalCounts { roi: &interval, name: "", cnts: RawCounts { nuc: &counts, coverage: 1 } }];
        let content = dummy.results();
        assert_eq!(content, expected);
    }
}
