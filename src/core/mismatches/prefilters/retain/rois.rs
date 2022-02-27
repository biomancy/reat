use std::borrow::Borrow;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::ops::Range;

use bio_types::genome::AbstractInterval;
use bio_types::genome::Interval;
use bio_types::genome::Position;

use crate::core::mismatches::prefilters::retain::ROIRetainer;
use crate::core::mismatches::roi::ROIMismatches;
use crate::core::workload::ROI;

// https://users.rust-lang.org/t/using-hashset-contains-with-tuple-types-without-takeing-ownership-of-the-values/65455/4
trait ROIHash {
    fn name(&self) -> &str;
    fn interval(&self) -> &Interval;
}

impl Hash for dyn ROIHash + '_ {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name().hash(state);
        self.interval().hash(state);
    }
}

impl PartialEq for dyn ROIHash + '_ {
    fn eq(&self, other: &Self) -> bool {
        self.name() == other.name() && self.interval() == other.interval()
    }
}
impl Eq for dyn ROIHash + '_ {}

impl ROIHash for (Interval, String) {
    fn name(&self) -> &str {
        &self.1
    }

    fn interval(&self) -> &Interval {
        &self.0
    }
}

impl ROIHash for (&Interval, &str) {
    fn name(&self) -> &str {
        self.1
    }

    fn interval(&self) -> &Interval {
        self.0
    }
}

impl<'a> Borrow<dyn ROIHash + 'a> for (Interval, String) {
    fn borrow(&self) -> &(dyn ROIHash + 'a) {
        self
    }
}

pub struct RetainROIFromList {
    hash: HashSet<(Interval, String)>,
}

impl RetainROIFromList {
    pub fn new(rois: Vec<ROI>) -> Self {
        let mut hash = HashSet::new();
        for r in rois.into_iter() {
            let (_, interval, name) = r.dissolve();
            hash.insert((interval, name));
        }
        Self { hash }
    }
}

impl ROIRetainer for RetainROIFromList {
    fn filter(&self, interval: &Interval, name: &str) -> bool {
        self.hash.contains::<dyn ROIHash>(&(interval, name))
    }
}
