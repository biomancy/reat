use crate::core::dna::NucCounts;
use crate::core::rpileup::ncounters::AggregatedNucCountsItem;

impl<'a, ItemInfo> AggregatedNucCountsItem<'a, ItemInfo> {
    pub fn seqnuc<'b>(&'a self, buffer: &'b mut Vec<NucCounts>) -> Option<&'a [NucCounts]> {
        // Gather sequenced nucleotides in each position
        match (self.forward, self.reverse, self.unstranded) {
            (Some(c), None, None) => Some(c),
            (None, Some(c), None) => Some(c),
            (None, None, Some(c)) => Some(c),
            (_, _, _) => {
                buffer.clear();
                buffer.resize(self.range.end as usize - self.range.start as usize, Default::default());
                for x in [self.forward, self.reverse, self.unstranded] {
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
