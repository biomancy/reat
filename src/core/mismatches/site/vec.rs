use std::io::Write;
use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use csv::Writer;
use itertools::{izip, Itertools};
use serde::ser::{SerializeSeq, SerializeStruct};
use serde::{Serialize, Serializer};

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::site::SiteDataRef;
use crate::core::mismatches::MismatchesVec;

use super::data::SiteDataVec;

#[derive(Clone)]
pub struct SiteMismatchesVec {
    contig: String,
    strand: Strand,
    pub data: SiteDataVec,
}

impl SiteMismatchesVec {
    pub fn new(contig: String, strand: Strand, data: SiteDataVec) -> Self {
        Self { contig, strand, data }
    }
}

impl MismatchesVec for SiteMismatchesVec {
    fn len(&self) -> usize {
        self.data.len()
    }

    fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    fn to_csv<F: Write>(&self, writer: &mut Writer<F>) -> csv::Result<()> {
        for data in &self.data {
            writer.serialize(&SerializeSiteRef { contig: &self.contig, strand: self.strand, data })?;
        }
        Ok(())
    }
}

struct SerializeSiteRef<'a> {
    contig: &'a str,
    strand: Strand,
    data: SiteDataRef<'a>,
}

impl Serialize for SerializeSiteRef<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("SiteMismatches", 9)?;
        state.serialize_field("contig", self.contig)?;
        state.serialize_field("pos", &self.data.pos)?;
        state.serialize_field("trstrand", self.strand.strand_symbol())?;
        state.serialize_field("refnuc", self.data.refnuc.symbol())?;
        state.serialize_field("prednuc", self.data.prednuc.symbol())?;
        state.serialize_field("A", &self.data.sequenced.A)?;
        state.serialize_field("C", &self.data.sequenced.C)?;
        state.serialize_field("G", &self.data.sequenced.G)?;
        state.serialize_field("T", &self.data.sequenced.T)?;
        state.end()
    }
}

#[cfg(test)]
mod test {
    use serde_test::{assert_ser_tokens, Token};

    use crate::core::mismatches::site::SiteData;

    use super::*;

    #[test]
    fn loci() {
        let data = SiteDataRef {
            pos: &13,
            refnuc: &Nucleotide::A,
            prednuc: &Nucleotide::G,
            sequenced: &NucCounts::new(1, 2, 3, 4),
        };
        assert_ser_tokens(
            &SerializeSiteRef { contig: "MySuperContig", strand: Strand::Unknown, data },
            &[
                Token::Struct { name: "SiteMismatches", len: 9 },
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
