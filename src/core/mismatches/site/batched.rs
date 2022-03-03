use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use itertools::{izip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::site::flat::REATSiteMismatches;
use crate::core::mismatches::BatchedMismatches;
use crate::core::mismatches::{utils, StrandingCounts};
use crate::core::strandutil::StrandedData;

use super::BatchedSiteMismatches;

#[derive(Clone)]
pub struct REATBatchedSiteMismatches {
    contig: String,
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

        let blocks = Self::infer_blocks(&loci);
        Self { contig, strand, loci, refnuc, prednuc, sequenced, blocks }
    }

    fn infer_blocks(loci: &[Position]) -> Vec<Range<Position>> {
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
        blocks
    }
}

impl AbstractInterval for REATBatchedSiteMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        *self.loci.first().unwrap()..self.loci.last().unwrap() + 1
    }
}

impl BatchedMismatches for REATBatchedSiteMismatches {
    type Flattened = REATSiteMismatches;

    fn trstrand(&self) -> Strand {
        self.strand
    }

    fn filter(mut self, mask: Vec<bool>) -> Self {
        debug_assert!([mask.len(), self.loci.len(), self.refnuc.len(), self.prednuc.len(), self.sequenced.len()]
            .iter()
            .all_equal());

        self.loci = utils::maskvec(self.loci, &mask);
        self.refnuc = utils::maskvec(self.refnuc, &mask);
        self.prednuc = utils::maskvec(self.prednuc, &mask);
        self.sequenced = utils::maskvec(self.sequenced, &mask);
        self.blocks = Self::infer_blocks(&self.loci);
        self
    }

    fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>> {
        let loci = utils::select_strands(self.loci, &strands, &cnts);
        let refnuc = utils::select_strands(self.refnuc, &strands, &cnts);
        let prednuc = utils::select_strands(self.prednuc, &strands, &cnts);
        let sequenced = utils::select_strands(self.sequenced, &strands, &cnts);

        let mut result = StrandedData { unknown: None, forward: None, reverse: None };
        if cnts.forward > 0 {
            result.forward =
                Some(Self::new(self.contig.clone(), Strand::Forward, loci.0, refnuc.0, prednuc.0, sequenced.0));
        }
        if cnts.reverse > 0 {
            result.reverse =
                Some(Self::new(self.contig.clone(), Strand::Reverse, loci.1, refnuc.1, prednuc.1, sequenced.1));
        }
        if cnts.unknown > 0 {
            result.unknown = Some(Self::new(self.contig, Strand::Unknown, loci.2, refnuc.2, prednuc.2, sequenced.2));
        }
        result
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

impl BatchedSiteMismatches for REATBatchedSiteMismatches {
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
