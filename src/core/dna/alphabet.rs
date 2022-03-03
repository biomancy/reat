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

impl Nucleotide {
    pub fn symbol(&self) -> &str {
        match self {
            Nucleotide::A => "A",
            Nucleotide::C => "C",
            Nucleotide::G => "G",
            Nucleotide::T => "T",
            Nucleotide::Unknown => "N",
        }
    }
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.symbol())
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

// impl From<Nucleotide> for ReqNucleotide {
//     fn from(nuc: Nucleotide) -> Self {
//         match nuc {
//             Nucleotide::A => ReqNucleotide::A,
//             Nucleotide::C => ReqNucleotide::C,
//             Nucleotide::G => ReqNucleotide::G,
//             Nucleotide::T => ReqNucleotide::T,
//             Nucleotide::Unknown => {
//                 panic!("Can't cast Nucleotide::Unknown to the ReqNucleotide")
//             }
//         }
//     }
// }
