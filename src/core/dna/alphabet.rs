use std::fmt::{Display, Formatter};

#[derive(Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Debug)]
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

impl Default for Nucleotide {
    fn default() -> Self {
        Nucleotide::Unknown
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

impl TryFrom<Nucleotide> for ReqNucleotide {
    type Error = ();

    fn try_from(nuc: Nucleotide) -> Result<Self, Self::Error> {
        match nuc {
            Nucleotide::A => Ok(ReqNucleotide::A),
            Nucleotide::C => Ok(ReqNucleotide::C),
            Nucleotide::G => Ok(ReqNucleotide::G),
            Nucleotide::T => Ok(ReqNucleotide::T),
            Nucleotide::Unknown => Err(()),
        }
    }
}

impl TryFrom<u8> for ReqNucleotide {
    type Error = ();

    fn try_from(nuc: u8) -> Result<Self, Self::Error> {
        match nuc {
            b'A' | b'a' => Ok(ReqNucleotide::A),
            b'C' | b'c' => Ok(ReqNucleotide::C),
            b'G' | b'g' => Ok(ReqNucleotide::G),
            b'T' | b't' => Ok(ReqNucleotide::T),
            _ => Err(()),
        }
    }
}
