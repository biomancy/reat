use bio_types::genome::Interval;
use bio_types::strand::Strand;

use crate::core::counting::NucCounts;
use crate::core::dna::Nucleotide;
use crate::core::summary::MismatchesSummary;

pub struct ROISummary {
    pub interval: Interval,
    pub strand: Strand,
    pub name: String,
    pub sequenced: NucCounts,
    pub mismatches: MismatchesSummary,
}

impl ROISummary {
    pub fn from_counts(
        interval: Interval,
        name: String,
        strand: Strand,
        sequence: &[Nucleotide],
        counts: &[NucCounts],
    ) -> Self {
        debug_assert_eq!(sequence.len(), counts.len());

        let mut mismatches = MismatchesSummary::zeros();
        let mut sequenced = NucCounts::zeros();
        for (seq, count) in sequence.iter().zip(counts) {
            match seq {
                Nucleotide::A => {
                    mismatches.A += *count;
                    sequenced.A += 1;
                }
                Nucleotide::C => {
                    mismatches.C += *count;
                    sequenced.C += 1;
                }
                Nucleotide::G => {
                    mismatches.G += *count;
                    sequenced.G += 1;
                }
                Nucleotide::T => {
                    mismatches.T += *count;
                    sequenced.T += 1;
                }
                _ => {}
            }
        }
        ROISummary { interval, strand, name, sequenced, mismatches }
    }
}

#[cfg(test)]
mod tests {
    use bio_types::strand::Same;

    use super::*;

    #[test]
    fn from_counts() {
        let interval = Interval::new("".into(), 0..5);
        let name = String::from("!");

        for strand in [Strand::Unknown, Strand::Reverse, Strand::Forward] {
            let summary = ROISummary::from_counts(
                interval.clone(),
                name.clone(),
                strand,
                &vec![Nucleotide::A, Nucleotide::Unknown, Nucleotide::C],
                &vec![
                    NucCounts { A: 10, C: 0, G: 15, T: 0 },
                    NucCounts { A: 10, C: 125, G: 15, T: 10 },
                    NucCounts { A: 9, C: 0, G: 8, T: 0 },
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
