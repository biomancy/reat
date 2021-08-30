use bio_types::genome::Interval;
use bio_types::strand::Strand;

use crate::core::counting::LocusCounts;
use crate::core::dna::Nucleotide;
use crate::core::summary::MismatchesSummary;

pub struct IntervalSummary {
    pub interval: Interval,
    pub strand: Strand,
    pub name: String,
    pub mismatches: MismatchesSummary,
}

impl IntervalSummary {
    pub fn from_counts(
        interval: Interval,
        name: String,
        strand: Strand,
        sequence: &[Nucleotide],
        counts: &[LocusCounts],
    ) -> Self {
        debug_assert_eq!(sequence.len(), counts.len());

        let mut mismatches = MismatchesSummary::zeros();
        for (seq, count) in sequence.iter().zip(counts) {
            match seq {
                Nucleotide::A => mismatches.A += *count,
                Nucleotide::C => mismatches.C += *count,
                Nucleotide::G => mismatches.G += *count,
                Nucleotide::T => mismatches.T += *count,
                _ => {}
            }
        }
        IntervalSummary { interval, strand, name, mismatches }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio_types::strand::Same;

    #[test]
    fn from_counts() {
        let interval = Interval::new("".into(), 0..5);
        let name = String::from("!");

        for strand in [Strand::Unknown, Strand::Reverse, Strand::Forward] {
            let summary = IntervalSummary::from_counts(
                interval.clone(),
                name.clone(),
                strand,
                &vec![Nucleotide::A, Nucleotide::Unknown, Nucleotide::C],
                &vec![
                    LocusCounts { A: 10, C: 0, G: 15, T: 0 },
                    LocusCounts { A: 10, C: 125, G: 15, T: 10 },
                    LocusCounts { A: 9, C: 0, G: 8, T: 0 },
                ],
            );
            assert_eq!(summary.name, name);
            assert!(summary.strand.same(&strand));

            let mut mismatches = MismatchesSummary::zeros();
            mismatches.A.A = 10;
            mismatches.A.G = 15;
            mismatches.C.A = 9;
            mismatches.C.G = 8;

            assert_eq!(summary.mismatches, mismatches);
        }
    }
}
