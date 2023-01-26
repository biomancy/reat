use std::ops::Range;

use bio_types::genome::Position;
use dyn_clone::DynClone;

pub use autoref::AutoRef;
pub use vcf::VCFCorrectedReference;

use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

mod autoref;
mod vcf;

#[derive(Clone, Copy, Debug)]
pub enum PredNucleotide {
    Homozygous(Nucleotide),
    Heterozygous((Nucleotide, Nucleotide)),
}

impl PredNucleotide {
    pub fn symbol(&self) -> &str {
        match self {
            PredNucleotide::Homozygous(nuc) => nuc.symbol(),
            PredNucleotide::Heterozygous((n1, n2)) => match (n1, n2) {
                (Nucleotide::A, Nucleotide::A) => "A/A",
                (Nucleotide::A, Nucleotide::C) => "A/C",
                (Nucleotide::A, Nucleotide::G) => "A/G",
                (Nucleotide::A, Nucleotide::T) => "A/T",

                (Nucleotide::C, Nucleotide::A) => "A/C",
                (Nucleotide::C, Nucleotide::C) => "C/C",
                (Nucleotide::C, Nucleotide::G) => "C/G",
                (Nucleotide::C, Nucleotide::T) => "C/T",

                (Nucleotide::G, Nucleotide::A) => "A/G",
                (Nucleotide::G, Nucleotide::C) => "C/G",
                (Nucleotide::G, Nucleotide::G) => "G/G",
                (Nucleotide::G, Nucleotide::T) => "G/T",

                (Nucleotide::T, Nucleotide::A) => "A/T",
                (Nucleotide::T, Nucleotide::C) => "C/T",
                (Nucleotide::T, Nucleotide::G) => "G/T",
                (Nucleotide::T, Nucleotide::T) => "T/T",
                _ => "N",
            },
        }
    }
}

impl Default for PredNucleotide {
    fn default() -> Self {
        PredNucleotide::Homozygous(Nucleotide::Unknown)
    }
}

pub struct RefEngineResult<'a> {
    pub predicted: &'a [PredNucleotide],
    pub reference: &'a [Nucleotide],
}

pub trait RefEngine: Send + DynClone {
    fn run(&mut self, contig: &str, range: Range<Position>, sequenced: &[NucCounts]);
    fn results(&self) -> RefEngineResult<'_>;
}
dyn_clone::clone_trait_object!(RefEngine);
