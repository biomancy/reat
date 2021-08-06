use bio_types::genome::{Interval, Locus};
use bio_types::strand::Strand;

use crate::rada::modules::counting::LocusCounts;
use crate::rada::modules::dna::ReqNucleotide;
use crate::rada::modules::summarization::MismatchSummary;

use super::StrandPredictor;
use crate::rada::modules::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};

pub struct StrandByAtoIEditing {
    min_mismatches: u32,
    min_freq: f32,
}

impl StrandByAtoIEditing {
    pub fn new(min_mismatches: u32, min_freq: f32) -> Self {
        Self {
            min_mismatches,
            min_freq,
        }
    }

    fn is_ok(&self, matches: u32, mismatches: u32) -> bool {
        let coverage = mismatches + matches;
        return coverage > 0
            && mismatches > self.min_mismatches
            && (mismatches as f32 / coverage as f32) > self.min_freq;
    }
}

impl StrandPredictor for StrandByAtoIEditing {}

impl IntervalStrandPredictor for StrandByAtoIEditing {
    fn predict(&self, _: &Interval, mismatches: &MismatchSummary) -> Strand {
        let a2g = self.is_ok(mismatches.A.A, mismatches.A.G);
        let t2c = self.is_ok(mismatches.T.T, mismatches.T.C);

        match (a2g, t2c) {
            (false, false) => Strand::Unknown,
            (false, true) => Strand::Reverse,
            (true, false) => Strand::Forward,
            (true, true) => {
                let a2g_coverage = mismatches.A.A + mismatches.A.G;
                let t2c_coverage = mismatches.T.T + mismatches.T.C;

                if a2g_coverage == 0 && t2c_coverage == 0 {
                    return Strand::Unknown;
                }

                let a2g = mismatches.A.G as f32 / a2g_coverage as f32;
                let t2c = mismatches.T.C as f32 / t2c_coverage as f32;
                if mismatches.A.G > 0 && a2g > t2c {
                    Strand::Forward
                } else if mismatches.T.C > 0 && t2c > a2g {
                    Strand::Reverse
                } else {
                    Strand::Unknown
                }
            }
        }
    }
}

impl LocusStrandPredictor for StrandByAtoIEditing {
    fn predict(&self, _: &Locus, refnuc: &ReqNucleotide, sequenced: &LocusCounts) -> Strand {
        match refnuc {
            ReqNucleotide::A => {
                if self.is_ok(sequenced.A, sequenced.G) {
                    Strand::Forward
                } else {
                    Strand::Reverse
                }
            }
            ReqNucleotide::C => Strand::Unknown,
            ReqNucleotide::G => Strand::Unknown,
            ReqNucleotide::T => {
                if self.is_ok(sequenced.T, sequenced.C) {
                    Strand::Reverse
                } else {
                    Strand::Forward
                }
            }
        }
    }
}
