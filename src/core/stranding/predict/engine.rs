use crate::core::mismatches::BatchedMismatches;
use crate::core::stranding::predict::StrandingContext;

use super::StrandingAlgo;
use super::StrandingEngine;

pub struct REATStrandingEngine<T: BatchedMismatches> {
    algo: Vec<Box<dyn StrandingAlgo<T>>>,
}

impl<T: BatchedMismatches> REATStrandingEngine<T> {
    pub fn new() -> Self {
        Self { algo: Vec::with_capacity(2) }
    }
}

impl<T: BatchedMismatches> StrandingEngine<T> for REATStrandingEngine<T> {
    fn add_algo(&mut self, algo: Box<dyn StrandingAlgo<T>>) {
        self.algo.push(algo)
    }

    fn strand(&self, items: Vec<T>) -> Vec<T> {
        let mut ctx = StrandingContext::new(items);

        for algo in &self.algo {
            ctx = algo.predict(ctx);
        }

        let (mut unstranded, stranded) = ctx.dissolve();
        unstranded.extend(stranded.into_iter());
        unstranded
    }
}
