use std::collections::HashMap;
use std::ops::Range;

pub use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Position};
use itertools::Itertools;

use super::SitesRetainer;

#[derive(Clone)]
pub struct RetainSitesFromIntervals {
    index: HashMap<String, IntervalTree<Position, ()>>,
}

impl RetainSitesFromIntervals {
    pub fn new(include: Vec<impl AbstractInterval>) -> Self {
        let mut index: HashMap<String, IntervalTree<Position, ()>> = HashMap::new();
        for record in include {
            if !index.contains_key(record.contig()) {
                index.insert(record.contig().into(), Default::default());
            }

            index.get_mut(record.contig()).map(|x| {
                // There must be no overlapping records!
                debug_assert!(x.find(record.range()).collect_vec().is_empty());
                x.insert(record.range(), ());
            });
        }
        Self { index }
    }
}

impl SitesRetainer for RetainSitesFromIntervals {
    #[inline]
    fn retained(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>> {
        match self.index.get(contig) {
            None => vec![],
            Some(map) => {
                let mut results = Vec::new();
                for hit in map.find(range) {
                    results.push(hit.interval().start..hit.interval().end)
                }
                results
            }
        }
    }
}
