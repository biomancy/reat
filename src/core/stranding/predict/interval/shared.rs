use bio_types::genome::{AbstractInterval, AbstractLocus, Interval};
use bio_types::strand::{Same, Strand};
use itertools::{zip, Itertools};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::interval::IntermediateIntervalMismatches;
use crate::core::mismatches::IntermediateMismatches;
use crate::core::stranding::predict::shared::StrandByAtoIEditing;
use crate::core::stranding::predict::shared::StrandByGenomicAnnotation;

use super::IntervalStrandPredictor;

impl StrandByAtoIEditing {
    #[inline]
    fn nuc_predict(&self, sequenced: &NucCounts, refnuc: &Nucleotide) -> Strand {
        match refnuc {
            Nucleotide::A => {
                if self.edited(sequenced.A, sequenced.G) {
                    Strand::Forward
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::T => {
                if self.edited(sequenced.T, sequenced.C) {
                    Strand::Reverse
                } else {
                    Strand::Unknown
                }
            }
            Nucleotide::C | Nucleotide::G | Nucleotide::Unknown => Strand::Unknown,
        }
    }
}

impl<T: IntermediateIntervalMismatches> IntervalStrandPredictor<T> for StrandByAtoIEditing {
    fn predict(&self, blocks: Vec<T>) -> Vec<T> {
        let mut results = Vec::with_capacity(blocks.len());

        let mut strands = Vec::new();
        let mut indices = Vec::new();
        for b in blocks {
            strands.clear();
            indices.clear();

            // Loop and safe strand for each site inside the interval
            for (sequenced, refnuc) in zip(b.ncounts(), b.refnuc()) {
                strands.push(self.nuc_predict(sequenced, refnuc));
            }
            // Split strands into constant regions
            let mut lastr = &strands[0];
            let mut cache = Vec::new();
            for (ind, strnd) in strands.iter().enumerate() {
                if !strnd.same(lastr) {
                    indices.push(ind);
                    cache.push(lastr);
                    lastr = strnd;
                }
            }
            cache.push(lastr);

            let mut b = b.split(&indices);
            debug_assert_eq!(b.len(), cache.len());
            for (mut x, strnd) in zip(b, cache) {
                x.set_strand(*strnd);
                results.push(x);
            }
        }
        results
    }
}

impl<T: IntermediateIntervalMismatches> IntervalStrandPredictor<T> for StrandByGenomicAnnotation {
    fn predict(&self, blocks: Vec<T>) -> Vec<T> {
        if blocks.is_empty() {
            return blocks;
        }

        let mut result = Vec::with_capacity(blocks.len());
        for mut i in blocks {
            debug_assert!(i.strand().is_unknown());

            // Split region into sub intervals with constant annotation if needed
            let interval = i.interval();
            let supfeatures = self.intervals_in(interval);
            debug_assert!(
                supfeatures.len() >= 1
                    && supfeatures.first().unwrap().range().start == i.interval().range().start
                    && supfeatures.last().unwrap().range().end == i.interval().range().end
            );

            if supfeatures.len() == 1 {
                let strand = self.predict(interval);
                i.set_strand(strand);
                result.push(i);
                continue;
            }

            let interst = interval.range().start;
            let splits = supfeatures.iter().skip(1).map(|x| (x.range().start - interst) as usize).collect_vec();
            let splited = i.split(&splits);

            debug_assert!(splited.iter().map(|x| x.interval()).zip(supfeatures.iter()).all(|(x, y)| x == y));

            for mut b in splited {
                b.set_strand(self.predict(b.interval()));
                result.push(b);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Neg;

    use bio_types::genome::{AbstractInterval, Interval};
    use bio_types::strand::Same;
    use bio_types::strand::Strand;
    use itertools::izip;

    use crate::core::dna::{NucCounts, Nucleotide};
    use crate::core::mismatches::interval::{BorrowedIntervalMismatches, IntervalMismatches, MockIntervalMismatches};
    use crate::core::stranding::predict::shared::StrandByAtoIEditing;

    use super::IntervalStrandPredictor;

    #[test]
    fn nuc_predict() {
        let dummy = StrandByAtoIEditing::new(10, 0.1);

        // Unknown strand
        for (sequenced, refnuc) in [
            (NucCounts::T(123), Nucleotide::C),
            (NucCounts::G(234), Nucleotide::Unknown),
            (NucCounts::A(32), Nucleotide::G),
            (NucCounts::new(170, 0, 170, 0), Nucleotide::C),
            (NucCounts::new(10, 0, 9, 0), Nucleotide::A),
            (NucCounts::new(0, 200, 0, 200), Nucleotide::G),
        ] {
            assert!(dummy.nuc_predict(&sequenced, &refnuc).is_unknown());
        }

        let dummy = StrandByAtoIEditing::new(8, 0.05);

        // Inferred strand
        for (matches, mismatches, strand) in
            [(8, 8, Strand::Forward), (100, 4, Strand::Unknown), (1, 7, Strand::Unknown), (10, 10, Strand::Forward)]
        {
            let cnts = NucCounts::new(matches, 0, mismatches, 0);
            assert!(dummy.nuc_predict(&cnts, &Nucleotide::A).same(&strand));

            let cnts = NucCounts::new(0, mismatches, 0, matches);
            assert!(dummy.nuc_predict(&cnts, &Nucleotide::T).same(&strand.neg()));
        }
    }

    #[test]
    fn intervals_by_editing() {
        let dummy = StrandByAtoIEditing::new(10, 0.5);
        let interval = Interval::new("".into(), 100..101);

        let assertok = |workload: Vec<(Vec<NucCounts>, Vec<Nucleotide>)>, strands, explen| {
            let mut blocks = Vec::new();
            for (sequenced, reference) in &workload {
                blocks.push(BorrowedIntervalMismatches::new(interval.clone(), Strand::Unknown, &reference, &sequenced))
            }

            let predicted = IntervalStrandPredictor::predict(&dummy, blocks);
            for (pred, expstrnd, explen) in izip!(predicted, strands, explen) {
                assert_eq!(pred.ncounts().len(), pred.refnuc().len());
                assert_eq!(pred.refnuc().len(), explen);
                assert!(pred.strand().same(&expstrnd));
            }
        };

        // Empty list
        let result: Vec<MockIntervalMismatches> = IntervalStrandPredictor::predict(&dummy, vec![]);
        assert!(result.is_empty());

        // Dummmies
        let a2g = |len| (vec![NucCounts { A: 5, C: 1, G: 20, T: 3 }].repeat(len), vec![Nucleotide::A].repeat(len));
        let t2c = |len| (vec![NucCounts { A: 2, C: 235, G: 4, T: 10 }].repeat(len), vec![Nucleotide::T].repeat(len));
        let unknown = |len| (vec![NucCounts { A: 1, C: 2, G: 3, T: 4 }].repeat(len), vec![Nucleotide::G].repeat(len));
        let chained = |items: Vec<(Vec<NucCounts>, Vec<Nucleotide>)>| {
            let mut result = (Vec::new(), Vec::new());
            for item in items {
                result.0.extend(item.0);
                result.1.extend(item.1);
            }
            result
        };

        // Single nucleotide
        assertok(
            vec![a2g(1), t2c(1), unknown(1)],
            vec![Strand::Forward, Strand::Reverse, Strand::Unknown],
            vec![1, 1, 1],
        );

        // Interval of nucleotides with each possible strand
        assertok(
            vec![a2g(10), t2c(3), unknown(150)],
            vec![Strand::Forward, Strand::Reverse, Strand::Unknown],
            vec![10, 3, 150],
        );

        // Intervals of nucleotides with multiple strands
        assertok(
            vec![
                chained(vec![a2g(13), t2c(49), a2g(1)]),
                chained(vec![a2g(1), t2c(12), unknown(199)]),
                chained(vec![t2c(123), unknown(1), a2g(1)]),
            ],
            vec![
                Strand::Forward,
                Strand::Reverse,
                Strand::Forward,
                Strand::Forward,
                Strand::Reverse,
                Strand::Unknown,
                Strand::Reverse,
                Strand::Unknown,
                Strand::Forward,
            ],
            vec![13, 49, 1, 1, 12, 199, 123, 1, 1],
        )
    }
}
