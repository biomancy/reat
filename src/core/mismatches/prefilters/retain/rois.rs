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
    fn contig(&self) -> &str;
    fn range(&self) -> &Range<Position>;
}

impl Hash for dyn ROIHash + '_ {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name().hash(state);
        self.range().hash(state);
        self.contig().hash(state);
    }
}

impl PartialEq for dyn ROIHash + '_ {
    fn eq(&self, other: &Self) -> bool {
        self.range() == other.range() && self.contig() == other.contig() && self.name() == other.name()
    }
}
impl Eq for dyn ROIHash + '_ {}

impl ROIHash for (String, Range<Position>, String) {
    fn name(&self) -> &str {
        &self.2
    }

    fn contig(&self) -> &str {
        &self.0
    }

    fn range(&self) -> &Range<Position> {
        &self.1
    }
}

impl ROIHash for (&str, &Range<Position>, &str) {
    fn name(&self) -> &str {
        self.2
    }

    fn contig(&self) -> &str {
        self.0
    }

    fn range(&self) -> &Range<Position> {
        self.1
    }
}

impl<'a> Borrow<dyn ROIHash + 'a> for (String, Range<Position>, String) {
    fn borrow(&self) -> &(dyn ROIHash + 'a) {
        self
    }
}

pub struct RetainROIFromList {
    hash: HashSet<(String, Range<Position>, String)>,
}

impl RetainROIFromList {
    pub fn new(rois: Vec<ROI>) -> Self {
        let mut hash = HashSet::new();
        for r in rois.into_iter() {
            let (_, interval, name) = r.dissolve();
            hash.insert((interval.contig().into(), interval.range(), name));
        }
        Self { hash }
    }
}

impl ROIRetainer for RetainROIFromList {
    #[inline]
    fn filter(&self, contig: &str, range: &Range<Position>, name: &str) -> bool {
        self.hash.contains::<dyn ROIHash>(&(contig, range, name))
    }
}
