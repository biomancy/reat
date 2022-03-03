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

// #[cfg(test)]
// mod test {
//     use bio_types::genome::Locus;
//     use bio_types::strand::Strand;
//
//     use crate::core::dna::NucCounts;
//     use crate::core::dna::Nucleotide;
//
//     use super::*;
//     use crate::core::mismatches::site::MockSiteMismatches;
//
//     #[test]
//     fn loci() {
//         let factory = |locus, strand, nuc, counts| {
//             let mut dummy = MockSiteMismatches::new();
//             dummy.expect_locus().return_const(locus);
//             dummy.expect_strand().return_const(strand);
//             dummy.expect_refnuc().return_const(nuc);
//             dummy.expect_ncounts().return_const(counts);
//             dummy
//         };
//
//         let workload = vec![
//             factory(Locus::new("1".into(), 10), Strand::Unknown, Nucleotide::Unknown, NucCounts::zeros()),
//             factory(Locus::new("2".into(), 20), Strand::Forward, Nucleotide::G, NucCounts::new(1, 2, 3, 4)),
//         ];
//         let mut saveto = Vec::new();
//         super::sites(&mut saveto, workload);
//
//         let result = String::from_utf8(saveto).unwrap();
//         let expected = "chr\tposition\tstrand\treference\tA\tC\tG\tT\n\
//                               1\t10\t.\tN\t0\t0\t0\t0\n\
//                               2\t20\t+\tG\t1\t2\t3\t4\n";
//         assert_eq!(&result, expected)
//     }
// }
