use bio_types::strand::Strand;
use derive_getters::Dissolve;

use crate::core::mismatches::BatchedMismatches;
use crate::core::stranding::predict::{StrandingAlgo, StrandingAlgoResult};

#[derive(Dissolve)]
pub struct StrandingContext<T: BatchedMismatches> {
    unstranded: Vec<T>,
    stranded: Vec<T>,
}

impl<T: BatchedMismatches> StrandingContext<T> {
    pub fn new(items: Vec<T>) -> Self {
        let (unstranded, stranded) = items.into_iter().partition(|x| x.trstrand().is_unknown());
        Self { unstranded, stranded }
    }

    pub fn apply(&mut self, algo: &Box<dyn StrandingAlgo<T>>) {
        self.unstranded = std::mem::take(&mut self.unstranded)
            .into_iter()
            .filter_map(|mut x| match algo.predict(&x) {
                StrandingAlgoResult::AllElements(strand) => match strand {
                    Strand::Forward | Strand::Reverse => {
                        x.restrand_all(strand);
                        self.stranded.push(x);
                        None
                    }
                    Strand::Unknown => Some(x),
                },
                StrandingAlgoResult::EachElement((strands, cnts)) => {
                    let x = x.restrand(strands, cnts);
                    x.forward.map(|e| self.stranded.push(e));
                    x.reverse.map(|e| self.stranded.push(e));
                    x.unknown
                }
            })
            .collect();
    }
}
