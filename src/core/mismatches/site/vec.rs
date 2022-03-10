use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use itertools::{izip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::site::flat::REATSiteMismatches;
use crate::core::mismatches::MismatchesVec;
use crate::core::mismatches::{utils, StrandingCounts};
use crate::core::strandutil::StrandedData;

use super::SiteMismatchesVec;

#[derive(Clone)]
pub struct REATSiteMismatchesVec {
    contig: String,
    strand: Strand,
    pos: Vec<Position>,
    refnuc: Vec<Nucleotide>,
    prednuc: Vec<Nucleotide>,
    sequenced: Vec<NucCounts>,
    blocks: Vec<Range<Position>>,
}

impl REATSiteMismatchesVec {
    pub fn new(
        contig: String,
        strand: Strand,
        pos: Vec<Position>,
        refnuc: Vec<Nucleotide>,
        prednuc: Vec<Nucleotide>,
        sequenced: Vec<NucCounts>,
    ) -> Self {
        // Not empty and equal size
        debug_assert!(!pos.is_empty());
        debug_assert!([pos.len(), refnuc.len(), prednuc.len(), sequenced.len()].iter().all_equal());
        // Sites must be ordered
        debug_assert!(pos.windows(2).all(|x| x[0] <= x[1]));

        let blocks = Self::infer_blocks(&pos);
        Self { contig, strand, pos, refnuc, prednuc, sequenced, blocks }
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

impl AbstractInterval for REATSiteMismatchesVec {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        *self.pos.first().unwrap()..self.pos.last().unwrap() + 1
    }
}

impl MismatchesVec for REATSiteMismatchesVec {
    type Flat = REATSiteMismatches;

    fn trstrand(&self) -> Strand {
        self.strand
    }

    // fn filter(mut self, mask: Vec<bool>) -> Self {
    //     debug_assert!([mask.len(), self.pos.len(), self.refnuc.len(), self.prednuc.len(), self.sequenced.len()]
    //         .iter()
    //         .all_equal());
    //
    //     self.pos = utils::maskvec(self.pos, &mask);
    //     self.refnuc = utils::maskvec(self.refnuc, &mask);
    //     self.prednuc = utils::maskvec(self.prednuc, &mask);
    //     self.sequenced = utils::maskvec(self.sequenced, &mask);
    //     self.blocks = Self::infer_blocks(&self.pos);
    //     self
    // }
    //
    // fn extend(&mut self, other: Self) {
    //     todo!()
    // }
    //
    // fn restrand(self, strands: Vec<Strand>, fwdto: &mut Self, revto: &mut Self) {
    //     todo!()
    // }

    // fn restrand(self, strands: Vec<Strand>, cnts: StrandingCounts) -> StrandedData<Option<Self>> {
    //     let mut pos = utils::select_strands(self.pos, &strands, &cnts);
    //     let mut refnuc = utils::select_strands(self.refnuc, &strands, &cnts);
    //     let mut prednuc = utils::select_strands(self.prednuc, &strands, &cnts);
    //     let mut sequenced = utils::select_strands(self.sequenced, &strands, &cnts);
    //
    //     let mut result = StrandedData { unknown: None, forward: None, reverse: None };
    //     for strand in [Strand::Forward, Strand::Reverse, Strand::Unknown] {
    //         if cnts[strand] > 0 {
    //             result[strand] = Some(Self::new(
    //                 self.contig.clone(),
    //                 strand,
    //                 std::mem::take(&mut pos[strand]),
    //                 std::mem::take(&mut refnuc[strand]),
    //                 std::mem::take(&mut prednuc[strand]),
    //                 std::mem::take(&mut sequenced[strand]),
    //             ));
    //         }
    //     }
    //     result
    // }

    // fn restrand_all(&mut self, strand: Strand) {
    //     self.strand = strand;
    // }

    fn flatten(self) -> Vec<Self::Flat> {
        izip!(self.pos.into_iter(), self.refnuc.into_iter(), self.prednuc.into_iter(), self.sequenced.into_iter())
            .map(|(pos, refn, predn, seq)| {
                REATSiteMismatches::new(self.contig.clone(), pos, self.strand, refn, predn, seq)
            })
            .collect()
    }
}

impl SiteMismatchesVec for REATSiteMismatchesVec {
    fn pos(&self) -> &[Position] {
        &self.pos
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
