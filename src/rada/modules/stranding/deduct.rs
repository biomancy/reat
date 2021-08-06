use bio_types::strand::ReqStrand;
use rust_htslib::bam::Record;

pub trait StrandDeductor {
    fn deduce(&self, record: &mut Record) -> ReqStrand;
}

pub enum StrandSpecificExperimentDesign {
    // Shortcuts (+ <- sequenced read orientation, - <- parent transcript orientation) => +-
    // Single end experiments
    Same, // forward read -> forward transcript (++, --)
    Flip, // forward read -> reverse transcript (+-, -+)
    // Paired end experiments
    Same1Flip2, // read1 same as the transcript, read2 flipped (++/--, +-/-+)
    Flip1Same2, // read1 flipped, read2 same as the transcript (+-/-+, ++/--)
}

pub struct DeductStrandByDesign {
    design: StrandSpecificExperimentDesign,
}

impl DeductStrandByDesign {
    pub fn new(design: StrandSpecificExperimentDesign) -> Self {
        DeductStrandByDesign { design }
    }

    fn flip(strand: &ReqStrand) -> ReqStrand {
        if *strand == ReqStrand::Forward {
            ReqStrand::Reverse
        } else {
            ReqStrand::Forward
        }
    }
}

impl StrandDeductor for DeductStrandByDesign {
    // MUTABLE IS NOT REALLY NEEDED HERE!
    fn deduce(&self, record: &mut Record) -> ReqStrand {
        let strand = record.strand();
        match self.design {
            StrandSpecificExperimentDesign::Same => strand,
            StrandSpecificExperimentDesign::Flip => DeductStrandByDesign::flip(&strand),
            StrandSpecificExperimentDesign::Same1Flip2 => {
                if record.is_first_in_template() {
                    strand
                } else {
                    DeductStrandByDesign::flip(&strand)
                }
            }
            StrandSpecificExperimentDesign::Flip1Same2 => {
                if record.is_first_in_template() {
                    DeductStrandByDesign::flip(&strand)
                } else {
                    strand
                }
            }
        }
    }
}