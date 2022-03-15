use std::ops::Range;

use bio_types::genome::Position;
use itertools::zip;

use crate::core::dna::NucCounts;
use crate::core::dna::{Nucleotide, ReqNucleotide};
use crate::core::io::fasta::FastaReader;
use crate::core::refpred::RefEngineResult;

use super::RefEngine;

#[derive(Clone)]
pub struct AutoRef<T: FastaReader> {
    mincoverage: u32,
    minfreq: f32,
    skip_hyperediting: bool,
    cache: Vec<Nucleotide>,
    reader: T,
}

impl<T: FastaReader> AutoRef<T> {
    pub fn new(mincoverage: u32, minfreq: f32, skip_hyperediting: bool, reader: T) -> Self {
        Self { mincoverage, minfreq, skip_hyperediting, cache: Vec::new(), reader }
    }

    #[inline]
    pub fn infer(&self, assembly: Nucleotide, sequenced: &NucCounts) -> Nucleotide {
        let coverage = sequenced.coverage();

        // if coverage is sufficient
        if coverage >= self.mincoverage {
            let (nuc, counts) = sequenced.mostfreq();
            // and the most abundant nucleotide is frequent enough
            if *counts as f32 / coverage as f32 >= self.minfreq {
                // check for a potential hyper-editing site
                let skip_hyperediting = self.skip_hyperediting
                    && ((assembly == Nucleotide::A && nuc == ReqNucleotide::G)
                        || (assembly == Nucleotide::T && nuc == ReqNucleotide::C));

                return if skip_hyperediting { assembly } else { nuc.into() };
            }
        }
        assembly
    }
}

impl<T: FastaReader> RefEngine for AutoRef<T> {
    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]) {
        self.cache.clear();
        self.cache.reserve(sequenced.len());

        self.reader.fetch(contig, range);
        let reference = self.reader.result();
        debug_assert!(reference.len() == sequenced.len());

        for (r, s) in zip(sequenced, reference) {
            let inferred = self.infer(*s, r);
            self.cache.push(inferred);
        }
    }

    fn results(&self) -> RefEngineResult<'_> {
        RefEngineResult { predicted: &self.cache, reference: &self.reader.result() }
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::{AbstractInterval, Interval};
    use mockall::Sequence;

    use crate::core::io::fasta::MockFastaReader;

    use super::*;

    #[test]
    fn infer() {
        let sequenced = NucCounts { A: 1, C: 2, G: 3, T: 4 };
        let assembly = Nucleotide::A;
        for (mincoverage, minfreq, result) in [
            (100, 0.0, assembly),
            (0, 1.0, assembly),
            (0, 0.41, assembly),
            (0, 0.0, Nucleotide::T),
            (0, 0.4, Nucleotide::T),
            (4, 0.4, Nucleotide::T),
        ] {
            let dummy = AutoRef::new(mincoverage, minfreq, false, MockFastaReader::new());
            assert_eq!(dummy.infer(assembly, &sequenced), result);
        }
    }

    #[test]
    fn results() {
        let intervals = vec![Interval::new("".into(), 1..4), Interval::new("chr1".into(), 100..105)];
        let sequenced = vec![
            (
                vec![NucCounts::A(1000), NucCounts::G(30), NucCounts::zeros()],
                vec![Nucleotide::G, Nucleotide::G, Nucleotide::A],
                vec![Nucleotide::A, Nucleotide::G, Nucleotide::A],
            ),
            (
                vec![NucCounts::G(12), NucCounts::G(10), NucCounts::C(32), NucCounts::T(16)],
                vec![Nucleotide::Unknown, Nucleotide::Unknown, Nucleotide::Unknown, Nucleotide::Unknown],
                vec![Nucleotide::G, Nucleotide::G, Nucleotide::C, Nucleotide::T],
            ),
        ];

        let mut reader = MockFastaReader::new();
        let mut seq = Sequence::new();
        for ind in 0..intervals.len() {
            reader.expect_fetch().once().return_const(()).in_sequence(&mut seq);
            reader.expect_result().once().return_const(sequenced[ind].1.clone()).in_sequence(&mut seq);
            reader.expect_result().once().return_const(sequenced[ind].1.clone()).in_sequence(&mut seq);
        }

        let mut dummy = AutoRef::new(10, 1f32, false, reader);

        for ind in 0..sequenced.len() {
            dummy.run(intervals[ind].contig(), intervals[ind].range(), &sequenced[ind].0);
            let result = dummy.results();
            assert_eq!(result.reference, sequenced[ind].1);
            assert_eq!(result.predicted, sequenced[ind].2);
        }
    }

    #[test]
    fn skip_hyper_editing() {
        let run = |expected, skip, sequenced, assembly| {
            for (ex, sk) in zip(expected, skip) {
                let dummy = AutoRef::new(0, 0f32, sk, MockFastaReader::new());
                assert_eq!(dummy.infer(assembly, sequenced), ex);
            }
        };

        let a2g = NucCounts { A: 1, C: 0, G: 99, T: 0 };
        run([Nucleotide::A, Nucleotide::G], [true, false], &a2g, Nucleotide::A);

        let t2c = NucCounts { A: 0, C: 3, G: 0, T: 1 };
        run([Nucleotide::T, Nucleotide::C], [true, false], &t2c, Nucleotide::T);
    }
}
