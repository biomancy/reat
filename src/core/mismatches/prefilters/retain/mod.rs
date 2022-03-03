use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;

pub use rois::RetainROIFromList;
pub use sites::RetainSitesFromIntervals;

mod rois;
mod sites;

pub trait SitesRetainer {
    fn filter(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>>;
}

pub trait ROIRetainer {
    fn filter(&self, contig: &str, range: &Range<Position>, strand: Strand, name: &str) -> bool;
}
