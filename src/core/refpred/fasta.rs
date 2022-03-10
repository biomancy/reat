use std::ops::Range;

use bio_types::genome::Position;

use crate::core::dna::NucCounts;
use crate::core::io::fasta::FastaReader;
use crate::core::refpred::RefEngine;

use super::RefEngineResult;

impl<T: FastaReader> RefEngine for T {
    fn run(&mut self, contig: &str, range: Range<Position>, _: &[NucCounts]) {
        self.fetch(contig, range)
    }

    fn results(&self) -> RefEngineResult<'_> {
        let result = self.result();
        RefEngineResult { predicted: result, reference: result }
    }
}
