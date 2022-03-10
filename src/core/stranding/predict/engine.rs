use bio_types::strand::Strand;

use crate::core::mismatches::MismatchesVec;
use crate::core::stranding::predict::StrandingAlgoResult;
use crate::core::strandutil::StrandedData;

use super::StrandingAlgo;
use super::StrandingEngine;

pub struct REATStrandingEngine<T> {
    algo: Vec<Box<dyn StrandingAlgo<T>>>,
}

impl<T> REATStrandingEngine<T> {
    pub fn new() -> Self {
        Self { algo: Vec::new() }
    }
    pub fn add(&mut self, algo: Box<dyn StrandingAlgo<T>>) {
        self.algo.push(algo)
    }
    pub fn clear(&mut self) {
        self.algo.clear()
    }
}

impl<T: MismatchesVec> StrandingEngine<T> for REATStrandingEngine<T> {
    fn strand(&self, mut items: StrandedData<Option<T>>) -> StrandedData<Option<T>> {
        if self.algo.is_empty() || items.unknown.is_none() {
            return items;
        }

        // for algo in &self.algo {
        //     items.unknown = items.unknown.and_then(|mut mm| match algo.predict(&mm) {
        //         // Identical strand for all elements
        //         StrandingAlgoResult::AllElements(strand) => match strand {
        //             // Pass through
        //             Strand::Unknown => Some(mm),
        //             // Restrand + extend
        //             Strand::Forward | Strand::Reverse => {
        //                 mm.restrand_all(strand);
        //
        //                 items[strand] = if let Some(mut other) = items[strand].take() {
        //                     other.extend(mm);
        //                     Some(other)
        //                 } else {
        //                     Some(mm)
        //                 };
        //
        //                 None
        //             }
        //         },
        //         StrandingAlgoResult::AllElements(Strand::Unknown) => Some(mm),
        //         StrandingAlgoResult::EachElement(_) => None,
        //     })

        // items.unknown = std::mem::take(&mut items.unknown).map(|x| {
        //     match algo.predict()
        // });
        //     .into_iter()
        //     .filter_map(|mut x| match algo.predict(&x) {
        //         StrandingAlgoResult::AllElements(strand) => match strand {
        //             Strand::Forward | Strand::Reverse => {
        //                 x.restrand_all(strand);
        //                 items[strand] = items[strand].map(||);
        //                 None
        //             }
        //             Strand::Unknown => Some(x),
        //         },
        //         StrandingAlgoResult::EachElement((strands, cnts)) => {
        //             let x = x.restrand(strands, cnts);
        //             x.forward.map(|e| items.forward.push(e));
        //             x.reverse.map(|e| items.reverse.push(e));
        //             x.unknown
        //         }
        //     })
        //     .collect();
        // }

        todo!()
    }
}

impl<T> Clone for REATStrandingEngine<T> {
    fn clone(&self) -> Self {
        Self { algo: self.algo.iter().map(|x| dyn_clone::clone_box(x.as_ref())).collect() }
    }
}
