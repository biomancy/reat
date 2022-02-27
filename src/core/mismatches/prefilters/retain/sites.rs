use std::collections::HashMap;
use std::ops::Range;

pub use bio::data_structures::interval_tree::IntervalTree;
use bio_types::genome::{AbstractInterval, Position};
use itertools::Itertools;

use super::SitesRetainer;

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

            let map = index.get_mut(record.contig()).unwrap();

            // There must be no overlapping records!
            debug_assert!(map.find(record.range()).collect_vec().is_empty());

            map.insert(record.range(), ());
        }
        Self { index }
    }
}

impl SitesRetainer for RetainSitesFromIntervals {
    fn filter(&self, contig: &str, range: Range<Position>) -> Vec<Range<Position>> {
        let map = self.index.get(contig).unwrap();

        let mut results = Vec::with_capacity(10);
        for hit in map.find(range) {
            results.push(hit.interval().start..hit.interval().end)
        }
        results
    }
}

// impl<T: IntermediateIntervalMismatches> AlwaysRetainFilter<T> for SitesFromIntervals {
//     fn partition(&self, objects: Vec<T>) -> (Vec<T>, Vec<T>) {
//         let (mut retain, mut tofilter) = (vec![], vec![]);
//
//         for obj in objects {
//             let map = self.index.get(obj.interval().contig()).unwrap();
//
//             let interstart = obj.interval().range().start;
//             let mut splits = vec![];
//             for hit in map.find(obj.interval().range()) {
//                 splits.push((hit.interval().start - interstart) as usize);
//                 splits.push((hit.interval().end - interstart) as usize);
//             }
//             if splits.len() == 0 {
//                 tofilter.push(obj);
//                 continue;
//             }
//
//             let splitted = obj.split(&splits);
//             assert!(splitted.len() >= 3);
//             for (ind, s) in splitted.into_iter().enumerate() {
//                 if ind % 2 == 0 {
//                     tofilter.push(s);
//                 } else {
//                     retain.push(s);
//                 }
//             }
//         }
//         (retain, tofilter)
//     }
// }
