use bio_types::strand::Strand;

use crate::core::mismatches::MismatchesVec;
use crate::core::strandutil::Stranded;

use super::StrandingAlgo;
use super::StrandingEngine;

#[derive(Default)]
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
    fn strand(&self, contig: &str, mut items: Stranded<T>) -> Stranded<T> {
        if self.algo.is_empty() || items.unknown.is_empty() {
            return items;
        }

        for algo in &self.algo {
            if items.unknown.is_empty() {
                return items;
            }
            algo.predict(contig, &mut items);
        }
        items
    }
}

impl<T> Clone for REATStrandingEngine<T> {
    fn clone(&self) -> Self {
        Self { algo: self.algo.iter().map(|x| dyn_clone::clone_box(x.as_ref())).collect() }
    }
}
