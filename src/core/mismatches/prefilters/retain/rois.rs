use std::borrow::Borrow;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::ops::Range;

use bio_types::genome::AbstractInterval;
use bio_types::genome::Position;
use bio_types::strand::Strand;

use crate::core::io::bed::BedRecord;
use crate::core::mismatches::prefilters::retain::ROIRetainer;

// https://users.rust-lang.org/t/using-hashset-contains-with-tuple-types-without-takeing-ownership-of-the-values/65455/4
trait ROIHash {
    fn contig(&self) -> &str;
    fn range(&self) -> &Range<Position>;
    fn strand(&self) -> &str;
    fn name(&self) -> &str;
}

impl Hash for dyn ROIHash + '_ {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name().hash(state);
        self.range().hash(state);
        self.contig().hash(state);
        self.strand().hash(state);
    }
}

impl PartialEq for dyn ROIHash + '_ {
    fn eq(&self, other: &Self) -> bool {
        self.range() == other.range()
            && self.strand() == other.strand()
            && self.contig() == other.contig()
            && self.name() == other.name()
    }
}

impl Eq for dyn ROIHash + '_ {}

impl ROIHash for (String, Range<Position>, String, String) {
    fn contig(&self) -> &str {
        &self.0
    }

    fn range(&self) -> &Range<Position> {
        &self.1
    }

    fn strand(&self) -> &str {
        &self.2
    }

    fn name(&self) -> &str {
        &self.3
    }
}

impl ROIHash for (&str, &Range<Position>, &str, &str) {
    fn contig(&self) -> &str {
        self.0
    }

    fn range(&self) -> &Range<Position> {
        self.1
    }

    fn strand(&self) -> &str {
        self.2
    }

    fn name(&self) -> &str {
        self.3
    }
}

impl<'a> Borrow<dyn ROIHash + 'a> for (String, Range<Position>, String, String) {
    fn borrow(&self) -> &(dyn ROIHash + 'a) {
        self
    }
}

#[derive(Clone)]
pub struct RetainROIFromList {
    hash: HashSet<(String, Range<Position>, String, String)>,
}

impl RetainROIFromList {
    pub fn new(rois: Vec<BedRecord>) -> Self {
        let mut hash = HashSet::new();
        for r in rois.into_iter() {
            let (name, strand, roi) = r.dissolve();
            hash.insert((roi.contig().into(), roi.range(), strand.strand_symbol().into(), name));
        }
        Self { hash }
    }
}

impl ROIRetainer for RetainROIFromList {
    #[inline]
    fn filter(&self, contig: &str, range: &Range<Position>, strand: Strand, name: &str) -> bool {
        self.hash.contains::<dyn ROIHash>(&(contig, range, strand.strand_symbol(), name))
    }
}
