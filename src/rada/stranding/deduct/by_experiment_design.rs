use std::fmt::{Display, Formatter};

use bio_types::strand::ReqStrand;
use derive_more::Constructor;

use crate::rada::read::AlignedRead;

use super::StrandDeductor;

#[derive(Copy, Clone, Eq, PartialEq)]
pub enum StrandSpecificExperimentDesign {
    // Notation (+ <- sequenced read orientation, - <- parent transcript orientation) => +-
    // Single end experiments
    Same, // forward read -> forward transcript (++, --)
    Flip, // forward read -> reverse transcript (+-, -+)
    // Paired end experiments
    Same1Flip2, // read1 same as the transcript, read2 flipped (++/--, +-/-+)
    Flip1Same2, // read1 flipped, read2 same as the transcript (+-/-+, ++/--)
}

impl Display for StrandSpecificExperimentDesign {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl StrandSpecificExperimentDesign {
    pub fn symbol(&self) -> &str {
        match self {
            StrandSpecificExperimentDesign::Same => "s",
            StrandSpecificExperimentDesign::Flip => "f",
            StrandSpecificExperimentDesign::Same1Flip2 => "sf",
            StrandSpecificExperimentDesign::Flip1Same2 => "fs",
        }
    }
}

#[derive(Constructor, Copy, Clone)]
pub struct DeductStrandByDesign {
    design: StrandSpecificExperimentDesign,
}

impl DeductStrandByDesign {
    fn flip(strand: &ReqStrand) -> ReqStrand {
        if *strand == ReqStrand::Forward {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        }
    }
}

impl<R: AlignedRead> StrandDeductor<R> for DeductStrandByDesign {
    fn deduce(&self, record: &R) -> ReqStrand {
        let strand = *record.strand();
        match self.design {
            StrandSpecificExperimentDesign::Same => strand,
            StrandSpecificExperimentDesign::Flip => DeductStrandByDesign::flip(&strand),
            StrandSpecificExperimentDesign::Same1Flip2 => {
                if record.is_first() {
                    strand
                } else {
                    DeductStrandByDesign::flip(&strand)
                }
            }
            StrandSpecificExperimentDesign::Flip1Same2 => {
                if record.is_first() {
                    DeductStrandByDesign::flip(&strand)
                } else {
                    strand
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::rada::read::MockRead;

    use super::*;

    #[test]
    fn flip() {
        assert_eq!(DeductStrandByDesign::flip(&ReqStrand::Forward), ReqStrand::Reverse);
        assert_eq!(DeductStrandByDesign::flip(&ReqStrand::Reverse), ReqStrand::Forward);
    }

    #[test]
    fn deduce_single_end() {
        let mut read = MockRead::new();
        for (design, strand, expected) in [
            (StrandSpecificExperimentDesign::Same, ReqStrand::Forward, ReqStrand::Forward),
            (StrandSpecificExperimentDesign::Same, ReqStrand::Reverse, ReqStrand::Reverse),
            (StrandSpecificExperimentDesign::Flip, ReqStrand::Forward, ReqStrand::Reverse),
            (StrandSpecificExperimentDesign::Flip, ReqStrand::Reverse, ReqStrand::Forward),
        ] {
            let dummy = DeductStrandByDesign::new(design);
            read.expect_strand().return_const(strand);
            assert_eq!(dummy.deduce(&read), expected);
            read.checkpoint();
        }
    }

    #[test]
    fn deduce_paired_end() {
        let mut read = MockRead::new();
        for (design, is_first, strand, expected) in [
            (StrandSpecificExperimentDesign::Same1Flip2, true, ReqStrand::Forward, ReqStrand::Forward),
            (StrandSpecificExperimentDesign::Same1Flip2, true, ReqStrand::Reverse, ReqStrand::Reverse),
            (StrandSpecificExperimentDesign::Same1Flip2, false, ReqStrand::Forward, ReqStrand::Reverse),
            (StrandSpecificExperimentDesign::Same1Flip2, false, ReqStrand::Reverse, ReqStrand::Forward),
            (StrandSpecificExperimentDesign::Flip1Same2, true, ReqStrand::Forward, ReqStrand::Reverse),
            (StrandSpecificExperimentDesign::Flip1Same2, true, ReqStrand::Reverse, ReqStrand::Forward),
            (StrandSpecificExperimentDesign::Flip1Same2, false, ReqStrand::Forward, ReqStrand::Forward),
            (StrandSpecificExperimentDesign::Flip1Same2, false, ReqStrand::Reverse, ReqStrand::Reverse),
        ] {
            let dummy = DeductStrandByDesign::new(design);
            read.expect_strand().return_const(strand);
            read.expect_is_first().return_const(is_first);
            assert_eq!(dummy.deduce(&read), expected);
            read.checkpoint();
        }
    }
}
