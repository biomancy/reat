use crate::rada::summarization::MismatchesSummary;
use bio_types::genome::Interval;
use bio_types::strand::Strand;

pub struct IntervalSummary {
    pub interval: Interval,
    pub strand: Strand,
    pub name: String,
    pub mismatches: MismatchesSummary,
}

impl IntervalSummary {
    pub fn zeros(interval: Interval, strand: Strand, name: String) -> Self {
        IntervalSummary { interval, strand, name, mismatches: MismatchesSummary::zeros() }
    }
}
