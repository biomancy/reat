use std::cmp::Ordering;
use std::io::Write;

use bio_types::strand::Strand;
use csv::Writer;
use itertools::Itertools;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::core::mismatches::roi::{ROIDataRef, ROIDataVec};
use crate::core::mismatches::MismatchesVec;

pub struct ROIMismatchesVec {
    contig: String,
    trstrand: Strand,
    pub data: ROIDataVec,
}

impl ROIMismatchesVec {
    pub fn new(contig: String, trstrand: Strand, data: ROIDataVec) -> Self {
        Self { contig, trstrand, data }
    }
}

impl MismatchesVec for ROIMismatchesVec {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn trstrand(&self) -> Strand {
        self.trstrand
    }

    fn len(&self) -> usize {
        self.data.len()
    }

    fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    fn ugly_in_contig_sort_and_to_csv<F: Write>(items: Vec<Self>, writer: &mut Writer<F>) -> csv::Result<()> {
        fn pos_then_strand_then_name(first: &SerializeROIRef, second: &SerializeROIRef) -> Ordering {
            let mut ord = first.data.roi.premasked.start.cmp(&second.data.roi.premasked.start);
            if ord.is_eq() {
                ord = first.data.roi.premasked.end.cmp(&second.data.roi.premasked.end);
            }
            if ord.is_eq() {
                ord = first.strand.strand_symbol().cmp(second.strand.strand_symbol());
            }
            if ord.is_eq() {
                ord = first.data.roi.name.cmp(second.data.roi.name);
            }
            ord
        }

        let iter = items
            .iter()
            .flat_map(|x| x.data.iter().map(|data| SerializeROIRef { contig: &x.contig, strand: x.trstrand, data }))
            .sorted_by(pos_then_strand_then_name);
        for item in iter {
            writer.serialize(item)?;
        }
        Ok(())
    }
}

struct SerializeROIRef<'a> {
    contig: &'a str,
    strand: Strand,
    data: ROIDataRef<'a>,
}

impl Serialize for SerializeROIRef<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("ROIMismatches", 29)?;
        state.serialize_field("contig", &self.contig)?;
        state.serialize_field("start", &self.data.roi.premasked.start)?;
        state.serialize_field("end", &self.data.roi.premasked.end)?;
        state.serialize_field("strand", &self.data.roi.strand.strand_symbol())?;
        state.serialize_field("name", &self.data.roi.name)?;
        state.serialize_field("trstrand", &self.strand.strand_symbol())?;
        state.serialize_field("coverage", &self.data.coverage)?;
        state.serialize_field("nucmasked", &self.data.roi.nucmasked())?;
        state.serialize_field("heterozygous", &self.data.heterozygous)?;
        state.serialize_field("#A", &self.data.prednuc.A)?;
        state.serialize_field("A->A", &self.data.mismatches.A.A)?;
        state.serialize_field("A->C", &self.data.mismatches.A.C)?;
        state.serialize_field("A->G", &self.data.mismatches.A.G)?;
        state.serialize_field("A->T", &self.data.mismatches.A.T)?;
        state.serialize_field("#C", &self.data.prednuc.C)?;
        state.serialize_field("C->A", &self.data.mismatches.C.A)?;
        state.serialize_field("C->C", &self.data.mismatches.C.C)?;
        state.serialize_field("C->G", &self.data.mismatches.C.G)?;
        state.serialize_field("C->T", &self.data.mismatches.C.T)?;
        state.serialize_field("#G", &self.data.prednuc.G)?;
        state.serialize_field("G->A", &self.data.mismatches.G.A)?;
        state.serialize_field("G->C", &self.data.mismatches.G.C)?;
        state.serialize_field("G->G", &self.data.mismatches.G.G)?;
        state.serialize_field("G->T", &self.data.mismatches.G.T)?;
        state.serialize_field("#T", &self.data.prednuc.T)?;
        state.serialize_field("T->A", &self.data.mismatches.T.A)?;
        state.serialize_field("T->C", &self.data.mismatches.T.C)?;
        state.serialize_field("T->G", &self.data.mismatches.T.G)?;
        state.serialize_field("T->T", &self.data.mismatches.T.T)?;
        state.end()
    }
}

#[cfg(test)]
mod test {
    use serde_test::{assert_ser_tokens, Token};

    use crate::core::dna::NucCounts;
    use crate::core::mismatches::roi::{NucMismatches, ROIDataRecordRef};

    use super::*;

    #[test]
    fn roi() {
        let record = ROIDataRecordRef {
            premasked: &(0..123),
            postmasked: &(1..100),
            subintervals: &vec![1..10, 20..100],
            name: &"MyRep".to_owned(),
            strand: &Strand::Forward,
        };
        let mm = NucMismatches {
            A: NucCounts::new(1, 2, 3, 4),
            C: NucCounts::new(5, 6, 7, 8),
            G: NucCounts::new(9, 10, 11, 12),
            T: NucCounts::new(13, 14, 15, 16),
        };
        let roi = ROIDataRef {
            roi: record,
            coverage: &13,
            prednuc: &NucCounts::new(1, 12, 3, 5),
            heterozygous: &13,
            mismatches: &mm,
        };

        assert_ser_tokens(
            &SerializeROIRef { contig: "chr1", strand: Strand::Unknown, data: roi },
            &[
                Token::Struct { name: "ROIMismatches", len: 29 },
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
                Token::U64(34),
                Token::Str("heterozygous"),
                Token::U64(13),
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
