use std::ops::Range;

use bio_types::annot::pos::Pos;
use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use itertools::{izip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::site::flat::REATSiteMismatches;
use crate::core::mismatches::{BatchedMismatches, Restranded};

use super::BinnedSiteMismatches;

#[derive(Clone)]
pub struct REATBatchedSiteMismatches {
    contig: String,
    range: Range<Position>,
    strand: Strand,
    loci: Vec<Position>,
    refnuc: Vec<Nucleotide>,
    prednuc: Vec<Nucleotide>,
    sequenced: Vec<NucCounts>,
    blocks: Vec<Range<Position>>,
}

impl REATBatchedSiteMismatches {
    pub fn new(
        contig: String,
        strand: Strand,
        loci: Vec<Position>,
        refnuc: Vec<Nucleotide>,
        prednuc: Vec<Nucleotide>,
        sequenced: Vec<NucCounts>,
    ) -> Self {
        // Not empty and equal size
        debug_assert!(!loci.is_empty());
        debug_assert!([loci.len(), refnuc.len(), prednuc.len(), sequenced.len()].iter().all_equal());
        // Sites must be ordered
        debug_assert!(loci.windows(2).all(|x| x[0] <= x[1]));

        let mut blocks = Vec::with_capacity(loci.len());
        let mut start = *loci.first().unwrap();
        let mut end = start + 1;
        for &next in loci[1..].iter() {
            if next == end {
                end += 1;
            } else {
                blocks.push(start..end);
                start = next;
                end = start + 1;
            }
        }
        blocks.push(start..end);

        let range = *loci.first().unwrap()..loci.last().unwrap() + 1;

        Self { contig, range, strand, loci, refnuc, prednuc, sequenced, blocks }
    }
}

impl AbstractInterval for REATBatchedSiteMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

impl BatchedMismatches for REATBatchedSiteMismatches {
    type Flattened = REATSiteMismatches;

    fn strand(&self) -> Strand {
        self.strand
    }

    fn filter(self, mask: Vec<bool>) -> Self {
        todo!()
    }

    fn restrand(self, strands: Vec<Strand>) -> Restranded<Self> {
        todo!()
    }

    fn restrand_all(&mut self, strand: Strand) {
        self.strand = strand;
    }

    fn flatten(self) -> Vec<Self::Flattened> {
        izip!(self.loci.into_iter(), self.refnuc.into_iter(), self.prednuc.into_iter(), self.sequenced.into_iter())
            .map(|(pos, refn, predn, seq)| {
                REATSiteMismatches::new(self.contig.clone(), pos, self.strand, refn, predn, seq)
            })
            .collect()
    }
}

impl BinnedSiteMismatches for REATBatchedSiteMismatches {
    fn pos(&self) -> &[Position] {
        &self.loci
    }

    fn refnuc(&self) -> &[Nucleotide] {
        &self.refnuc
    }

    fn prednuc(&self) -> &[Nucleotide] {
        &self.prednuc
    }

    fn seqnuc(&self) -> &[NucCounts] {
        &self.sequenced
    }

    fn blocks(&self) -> &[Range<Position>] {
        &self.blocks
    }
}
