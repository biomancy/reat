use bio_types::strand::Strand;
use derive_getters::Getters;
use derive_more::Constructor;

use crate::core::dna::Nucleotide;

#[derive(Constructor, Getters, Copy, Clone)]
pub struct StrandByAtoIEditing {
    minmismatches: u32,
    minfreq: f32,
}

impl StrandByAtoIEditing {
    #[inline]
    pub fn edited(&self, matches: u32, mismatches: u32) -> bool {
        let coverage = mismatches + matches;
        coverage > 0 && mismatches >= self.minmismatches && (mismatches as f32 / coverage as f32) >= self.minfreq
    }
}
