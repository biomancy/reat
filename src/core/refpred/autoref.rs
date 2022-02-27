use std::ops::Range;
use std::path::PathBuf;

use bio_types::genome::{Interval, Position};
use itertools::zip;
use rust_htslib::faidx;

use crate::core::dna::NucCounts;
use crate::core::dna::{Nucleotide, ReqNucleotide};
use crate::core::io::fasta::FastaReader;

use super::RefEngine;

pub struct AutoRef<T: FastaReader> {
    mincoverage: u32,
    minfreq: f32,
    skip_hyperediting: bool,
    flatcache: Vec<Nucleotide>,
    ranges: Vec<(usize, usize)>,
    reader: T,
}

impl<T: FastaReader> AutoRef<T> {
    pub fn new(mincoverage: u32, minfreq: f32, skip_hyperediting: bool, reader: T) -> Self {
        Self { mincoverage, minfreq, skip_hyperediting, flatcache: Vec::new(), ranges: Vec::new(), reader }
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
    fn reset(&mut self) {
        self.flatcache.clear();
        self.ranges.clear();
    }

    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]) {
        self.flatcache.reserve(self.flatcache.len() + sequenced.len());

        self.reader.fetch(contig, range);
        let reference = self.reader.result();

        let begin = self.flatcache.len();
        for (r, s) in zip(sequenced, reference) {
            let inferred = self.infer(*s, r);
            self.flatcache.push(inferred);
        }
        let end = self.flatcache.len();

        self.ranges.push((begin, end));
    }

    fn results(&self) -> &[Nucleotide] {
        self.ranges.iter().map(|(start, end)| &self.flatcache[*start..*end]).collect()
    }
}

impl<T: FastaReader + Clone> Clone for AutoRef<T> {
    fn clone(&self) -> Self {
        Self {
            mincoverage: self.mincoverage,
            minfreq: self.minfreq,
            skip_hyperediting: self.skip_hyperediting,
            flatcache: Vec::new(),
            ranges: Vec::with_capacity(100),
            reader: self.reader.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use bio_types::genome::AbstractInterval;
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
                vec![Nucleotide::Unknown, Nucleotide::Unknown, Nucleotide::Unknown],
                vec![Nucleotide::G, Nucleotide::G, Nucleotide::C],
            ),
        ];

        let mut reader = MockFastaReader::new();
        let mut seq = Sequence::new();
        for ind in 0..intervals.len() {
            reader.expect_fetch().once().return_const(()).in_sequence(&mut seq);
            reader.expect_result().once().return_const(sequenced[ind].1.clone()).in_sequence(&mut seq);
        }

        let mut dummy = AutoRef::new(10, 1f32, false, reader);

        dummy.reset();
        for ind in 0..sequenced.len() {
            dummy.run(intervals[ind].contig(), intervals[ind].range(), &sequenced[ind].0);
        }
        let result = dummy.results();
        assert_eq!(result.len(), sequenced.len());
        for ind in 0..sequenced.len() {
            assert_eq!(result[ind], sequenced[ind].2)
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
