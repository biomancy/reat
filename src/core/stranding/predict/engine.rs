use crate::core::mismatches::BatchedMismatches;
use crate::core::stranding::predict::StrandingContext;

use super::StrandingAlgo;
use super::StrandingEngine;

pub struct REATStrandingEngine<T> {
    algo: Vec<Box<dyn StrandingAlgo<T>>>,
}

impl<T> REATStrandingEngine<T> {
    pub fn new() -> Self {
        Self { algo: Vec::with_capacity(2) }
    }
    pub fn add_algo(&mut self, algo: Box<dyn StrandingAlgo<T>>) {
        self.algo.push(algo)
    }
    pub fn clear(&mut self) {
        self.algo.clear()
    }
}

impl<T: BatchedMismatches> StrandingEngine<T> for REATStrandingEngine<T> {
    fn strand(&self, items: Vec<T>) -> Vec<T> {
        if self.algo.is_empty() {
            return items;
        }

        let mut ctx = StrandingContext::new(items);

        for algo in &self.algo {
            ctx.apply(algo);
        }

        let (mut unstranded, stranded) = ctx.dissolve();
        unstranded.extend(stranded.into_iter());
        unstranded
    }
}

impl<T> Clone for REATStrandingEngine<T> {
    fn clone(&self) -> Self {
        Self { algo: self.algo.iter().map(|x| dyn_clone::clone_box(x.as_ref())).collect() }
    }
}
