use bio_types::genome::{AbstractLocus, Position};
use bio_types::strand::Strand;
use derive_more::Constructor;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::core::dna::{NucCounts, Nucleotide};

pub use super::SiteMismatches;

#[derive(Clone, Debug, Constructor)]
pub struct REATSiteMismatches {
    contig: String,
    pos: Position,
    strand: Strand,
    refnuc: Nucleotide,
    prednuc: Nucleotide,
    sequenced: NucCounts,
}

impl AbstractLocus for REATSiteMismatches {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn pos(&self) -> Position {
        self.pos
    }
}

impl SiteMismatches for REATSiteMismatches {
    fn trstrand(&self) -> Strand {
        self.strand
    }

    fn refnuc(&self) -> Nucleotide {
        self.refnuc
    }

    fn prednuc(&self) -> Nucleotide {
        self.prednuc
    }

    fn sequenced(&self) -> NucCounts {
        self.sequenced
    }
}

impl Serialize for REATSiteMismatches {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("REATSiteMismatches", 9)?;
        state.serialize_field("contig", self.contig())?;
        state.serialize_field("pos", &self.pos())?;
        state.serialize_field("trstrand", self.trstrand().strand_symbol())?;
        state.serialize_field("refnuc", self.refnuc.symbol())?;
        state.serialize_field("prednuc", self.prednuc.symbol())?;
        state.serialize_field("A", &self.sequenced.A)?;
        state.serialize_field("C", &self.sequenced.C)?;
        state.serialize_field("G", &self.sequenced.G)?;
        state.serialize_field("T", &self.sequenced.T)?;
        state.end()
    }
}

#[cfg(test)]
mod test {
    use serde_test::{assert_ser_tokens, Token};

    use super::*;

    #[test]
    fn loci() {
        let site = REATSiteMismatches {
            contig: "MySuperContig".to_string(),
            pos: 13,
            strand: Strand::Unknown,
            refnuc: Nucleotide::A,
            prednuc: Nucleotide::G,
            sequenced: NucCounts::new(1, 2, 3, 4),
        };
        assert_ser_tokens(
            &site,
            &[
                Token::Struct { name: "REATSiteMismatches", len: 9 },
                Token::Str("contig"),
                Token::Str("MySuperContig"),
                Token::Str("pos"),
                Token::U64(13),
                Token::Str("trstrand"),
                Token::Str("."),
                Token::Str("refnuc"),
                Token::Str("A"),
                Token::Str("prednuc"),
                Token::Str("G"),
                Token::Str("A"),
                Token::U32(1),
                Token::Str("C"),
                Token::U32(2),
                Token::Str("G"),
                Token::U32(3),
                Token::Str("T"),
                Token::U32(4),
                Token::StructEnd,
            ],
        );
    }
}
