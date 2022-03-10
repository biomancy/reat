use bio_types::genome::AbstractInterval;
use bio_types::strand::Strand;
use derive_more::Constructor;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::core::dna::NucCounts;
use crate::core::mismatches::roi::NucMismatches;
use crate::core::workload::ROI;

pub use super::ROIMismatches;

#[derive(Clone, Debug, Constructor)]
pub struct REATROIMismatches {
    roi: ROI,
    trstrand: Strand,
    coverage: u32,
    prednuc: NucCounts,
    mismatches: NucMismatches,
}

impl ROIMismatches for REATROIMismatches {
    fn roi(&self) -> &ROI {
        &self.roi
    }

    fn trstrand(&self) -> Strand {
        self.trstrand
    }

    fn coverage(&self) -> u32 {
        self.coverage
    }

    fn prednuc(&self) -> NucCounts {
        self.prednuc
    }

    fn mismatches(&self) -> &NucMismatches {
        &self.mismatches
    }
}

impl Serialize for REATROIMismatches {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("REATROIMismatches", 28)?;
        state.serialize_field("contig", self.roi.contig())?;
        state.serialize_field("start", &self.roi.premasked().start)?;
        state.serialize_field("end", &self.roi.premasked().end)?;
        state.serialize_field("strand", &self.roi.strand().strand_symbol())?;
        state.serialize_field("name", &self.roi.name())?;
        state.serialize_field("trstrand", &self.trstrand.strand_symbol())?;
        state.serialize_field("coverage", &self.coverage())?;
        state.serialize_field("nucmasked", &self.roi.nucmasked())?;
        state.serialize_field("#A", &self.prednuc().A)?;
        state.serialize_field("A->A", &self.mismatches().A.A)?;
        state.serialize_field("A->C", &self.mismatches().A.C)?;
        state.serialize_field("A->G", &self.mismatches().A.G)?;
        state.serialize_field("A->T", &self.mismatches().A.T)?;
        state.serialize_field("#C", &self.prednuc().C)?;
        state.serialize_field("C->A", &self.mismatches().C.A)?;
        state.serialize_field("C->C", &self.mismatches().C.C)?;
        state.serialize_field("C->G", &self.mismatches().C.G)?;
        state.serialize_field("C->T", &self.mismatches().C.T)?;
        state.serialize_field("#G", &self.prednuc().G)?;
        state.serialize_field("G->A", &self.mismatches().G.A)?;
        state.serialize_field("G->C", &self.mismatches().G.C)?;
        state.serialize_field("G->G", &self.mismatches().G.G)?;
        state.serialize_field("G->T", &self.mismatches().G.T)?;
        state.serialize_field("#T", &self.prednuc().T)?;
        state.serialize_field("T->A", &self.mismatches().T.A)?;
        state.serialize_field("T->C", &self.mismatches().T.C)?;
        state.serialize_field("T->G", &self.mismatches().T.G)?;
        state.serialize_field("T->T", &self.mismatches().T.T)?;
        state.end()
    }
}

#[cfg(test)]
mod test {
    use serde_test::{assert_ser_tokens, Token};

    use super::*;

    #[test]
    fn roi() {
        let ROI = ROI::new("chr1".into(), 0..123, vec![0..10, 100..123], "MyRep".into(), Strand::Forward);
        let mut mm = NucMismatches::zeros();
        mm.A = NucCounts::new(1, 2, 3, 4);
        mm.C = NucCounts::new(5, 6, 7, 8);
        mm.G = NucCounts::new(9, 10, 11, 12);
        mm.T = NucCounts::new(13, 14, 15, 16);
        let ROI = REATROIMismatches {
            roi: ROI,
            trstrand: Strand::Unknown,
            coverage: 13,
            prednuc: NucCounts::new(1, 12, 3, 5),
            mismatches: mm,
        };

        assert_ser_tokens(
            &ROI,
            &[
                Token::Struct { name: "REATROIMismatches", len: 28 },
                Token::Str("contig"),
                Token::Str("chr1"),
                Token::Str("start"),
                Token::U64(0),
                Token::Str("end"),
                Token::U64(123),
                Token::Str("strand"),
                Token::Str("+"),
                Token::Str("name"),
                Token::Str("MyRep"),
                Token::Str("trstrand"),
                Token::Str("."),
                Token::Str("coverage"),
                Token::U32(13),
                Token::Str("nucmasked"),
                Token::U64(90),
                Token::Str("#A"),
                Token::U32(1),
                Token::Str("A->A"),
                Token::U32(1),
                Token::Str("A->C"),
                Token::U32(2),
                Token::Str("A->G"),
                Token::U32(3),
                Token::Str("A->T"),
                Token::U32(4),
                Token::Str("#C"),
                Token::U32(12),
                Token::Str("C->A"),
                Token::U32(5),
                Token::Str("C->C"),
                Token::U32(6),
                Token::Str("C->G"),
                Token::U32(7),
                Token::Str("C->T"),
                Token::U32(8),
                Token::Str("#G"),
                Token::U32(3),
                Token::Str("G->A"),
                Token::U32(9),
                Token::Str("G->C"),
                Token::U32(10),
                Token::Str("G->G"),
                Token::U32(11),
                Token::Str("G->T"),
                Token::U32(12),
                Token::Str("#T"),
                Token::U32(5),
                Token::Str("T->A"),
                Token::U32(13),
                Token::Str("T->C"),
                Token::U32(14),
                Token::Str("T->G"),
                Token::U32(15),
                Token::Str("T->T"),
                Token::U32(16),
                Token::StructEnd,
            ],
        );
    }
}
