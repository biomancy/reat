use std::fmt::{Display, Formatter};

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
#[allow(non_snake_case)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    Unknown,
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Nucleotide::A => {
                write!(f, "A")
            }
            Nucleotide::C => {
                write!(f, "C")
            }
            Nucleotide::G => {
                write!(f, "G")
            }
            Nucleotide::T => {
                write!(f, "T")
            }
            Nucleotide::Unknown => {
                write!(f, "N")
            }
        }
    }
}

impl From<u8> for Nucleotide {
    fn from(symbol: u8) -> Self {
        match symbol {
            b'A' | b'a' => Nucleotide::A,
            b'C' | b'c' => Nucleotide::C,
            b'G' | b'g' => Nucleotide::G,
            b'T' | b't' => Nucleotide::T,
            _ => Nucleotide::Unknown,
        }
    }
}

impl From<ReqNucleotide> for Nucleotide {
    fn from(nuc: ReqNucleotide) -> Self {
        match nuc {
            ReqNucleotide::A => Nucleotide::A,
            ReqNucleotide::C => Nucleotide::C,
            ReqNucleotide::G => Nucleotide::G,
            ReqNucleotide::T => Nucleotide::T,
        }
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
#[allow(non_snake_case)]
pub enum ReqNucleotide {
    A,
    C,
    G,
    T,
}

impl From<Nucleotide> for ReqNucleotide {
    fn from(nuc: Nucleotide) -> Self {
        match nuc {
            Nucleotide::A => ReqNucleotide::A,
            Nucleotide::C => ReqNucleotide::C,
            Nucleotide::G => ReqNucleotide::G,
            Nucleotide::T => ReqNucleotide::T,
            Nucleotide::Unknown => {
                panic!("Can't cast Nucleotide::Unknown to the ReqNucleotide")
            }
        }
    }
}
