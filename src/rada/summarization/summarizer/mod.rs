use bio_types::genome::{AbstractInterval, Interval, Locus};
use bio_types::strand::Strand;
use itertools::izip;
use itertools::Itertools;
use rust_htslib::faidx;

pub use base::BaseSummarizer;

use crate::rada::counting::{CountsBufferContent, LocusCounts, NucCounterContent};
use crate::rada::dna::{Nucleotide, ReqNucleotide};
use crate::rada::refnuc::RefNucPredictor;
use crate::rada::stranding::predict::{IntervalStrandPredictor, LocusStrandPredictor};
use crate::rada::summarization::summary::{IntervalSummary, LocusSummary, MismatchesSummary};

mod base;

pub trait IntervalSummarizer {
    fn summarize(&self, name: &str, content: NucCounterContent) -> Vec<IntervalSummary>;
}

pub trait LocusSummarizer {
    fn summarize(&self, content: NucCounterContent) -> Vec<LocusSummary>;
}

pub trait Summarizer: IntervalSummarizer + LocusSummarizer {}
