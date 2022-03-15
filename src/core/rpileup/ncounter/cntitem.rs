use crate::core::dna::NucCounts;
use crate::core::rpileup::ncounter::{InnerNucCounts};

impl<'a, Data> InnerNucCounts<'a, Data> {
    pub fn seqnuc<'b>(&'a self, buffer: &'b mut Vec<NucCounts>) -> Option<&'a [NucCounts]> {
        // Gather sequenced nucleotides in each position
        match (self.cnts.forward, self.cnts.reverse, self.cnts.unknown) {
            (Some(c), None, None) => Some(c),
            (None, Some(c), None) => Some(c),
            (None, None, Some(c)) => Some(c),
            _ => {
                buffer.clear();
                buffer.resize(self.range.end as usize - self.range.start as usize, Default::default());
                for x in [self.cnts.forward, self.cnts.reverse, self.cnts.unknown] {
                    if let Some(x) = x {
                        debug_assert!(x.len() == buffer.len());
                        for (bc, xc) in buffer.iter_mut().zip(x.iter()) {
                            *bc += *xc;
                        }
                    }
                }
                None
            }
        }
    }
}
