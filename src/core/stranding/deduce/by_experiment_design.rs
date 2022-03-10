use std::fmt::{Display, Formatter};

use bio_types::strand::ReqStrand;
use derive_more::Constructor;

use crate::core::read::AlignedRead;

use super::StrandDeducer;

#[derive(Copy, Clone, Eq, PartialEq)]
pub enum StrandSpecificExperimentDesign {
    // Notation: aligned read strand/transcript strand
    // Single end experiments
    Same, // read strand = transcript strand -> (++, --)
    Flip, // read strand = reverse transcript strand -> (+-, -+)
    // Paired end experiments
    Same1Flip2, // read1 strand = transcript strand, read2 strand =  transcript strand -> (++/--, +-/-+)
    Flip1Same2, // read1 strand = reversed transcript strand, read2 strand = transcript strand -> (+-/-+, ++/--)
}

impl Display for StrandSpecificExperimentDesign {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl StrandSpecificExperimentDesign {
    pub fn symbol(&self) -> &str {
        match self {
            StrandSpecificExperimentDesign::Same => "s (++, --)",
            StrandSpecificExperimentDesign::Flip => "f (+-, -+)",
            StrandSpecificExperimentDesign::Same1Flip2 => "s/f (++/--, +-/-+)",
            StrandSpecificExperimentDesign::Flip1Same2 => "f/s (+-/-+, ++/--)",
        }
    }
}

#[derive(Constructor, Copy, Clone)]
pub struct DeduceStrandByDesign {
    design: StrandSpecificExperimentDesign,
}

impl DeduceStrandByDesign {
    fn flip(strand: &ReqStrand) -> ReqStrand {
        if *strand == ReqStrand::Forward {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        }
    }
}

impl<R: AlignedRead> StrandDeducer<R> for DeduceStrandByDesign {
    fn deduce(&self, record: &R) -> ReqStrand {
        let strand = *record.strand();
        match self.design {
            StrandSpecificExperimentDesign::Same => strand,
            StrandSpecificExperimentDesign::Flip => DeduceStrandByDesign::flip(&strand),
            StrandSpecificExperimentDesign::Same1Flip2 => {
                if record.is_first() {
                    strand
                } else {
                    DeduceStrandByDesign::flip(&strand)
                }
            }
            StrandSpecificExperimentDesign::Flip1Same2 => {
                if record.is_first() {
                    DeduceStrandByDesign::flip(&strand)
                } else {
                    strand
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::read::MockRead;

    use super::*;

    #[test]
    fn symbol_and_format() {
        use StrandSpecificExperimentDesign::*;
        for (symbol, design) in [
            ("s (++, --)", Same),
            ("f (+-, -+)", Flip),
            ("s/f (++/--, +-/-+)", Same1Flip2),
            ("f/s (+-/-+, ++/--)", Flip1Same2),
        ] {
            assert_eq!(design.symbol(), symbol);
            assert_eq!(format!("{}", design), symbol);
        }
    }

    #[test]
    fn flip() {
        assert_eq!(DeduceStrandByDesign::flip(&ReqStrand::Forward), ReqStrand::Reverse);
        assert_eq!(DeduceStrandByDesign::flip(&ReqStrand::Reverse), ReqStrand::Forward);
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
            let dummy = DeduceStrandByDesign::new(design);
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
            let dummy = DeduceStrandByDesign::new(design);
            read.expect_strand().return_const(strand);
            read.expect_is_first().return_const(is_first);
            assert_eq!(dummy.deduce(&read), expected);
            read.checkpoint();
        }
    }
}
