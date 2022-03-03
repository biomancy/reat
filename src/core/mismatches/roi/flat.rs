use bio_types::genome::AbstractInterval;
use bio_types::strand::Strand;
use derive_more::Constructor;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use crate::core::dna::NucCounts;
use crate::core::mismatches::roi::MismatchesSummary;
use crate::core::workload::ROI;

pub use super::ROIMismatches;

#[derive(Clone, Debug, Constructor)]
pub struct REATROIMismatches {
    roi: ROI,
    trstrand: Strand,
    coverage: u32,
    prednuc: NucCounts,
    mismatches: MismatchesSummary,
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

    fn mismatches(&self) -> &MismatchesSummary {
        &self.mismatches
    }
}

impl Serialize for REATROIMismatches {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("REATROIMismatches", 9)?;

        state.serialize_field("contig", self.roi.contig())?;
        state.serialize_field("start", &self.roi.interval().range().start)?;
        state.serialize_field("end", &self.roi.interval().range().end)?;
        state.serialize_field("strand", &self.roi.strand().strand_symbol())?;
        state.serialize_field("name", &self.roi.name())?;
        state.serialize_field("trstrand", &self.trstrand.strand_symbol())?;
        state.serialize_field("coverage", &self.coverage())?;
        state.serialize_field("nucmasked", &self.roi.nucmasked())?;
        state.serialize_field("#A", &self.prednuc().A)?;
        state.serialize_field("A->A", &self.mismatches().A.A)?;
        state.serialize_field("A->T", &self.mismatches().A.T)?;
        state.serialize_field("A->G", &self.mismatches().A.G)?;
        state.serialize_field("A->C", &self.mismatches().A.C)?;
        state.serialize_field("#T", &self.prednuc().T)?;
        state.serialize_field("T->A", &self.mismatches().T.A)?;
        state.serialize_field("T->T", &self.mismatches().T.T)?;
        state.serialize_field("T->G", &self.mismatches().T.G)?;
        state.serialize_field("T->C", &self.mismatches().T.C)?;
        state.serialize_field("#G", &self.prednuc().G)?;
        state.serialize_field("G->A", &self.mismatches().G.A)?;
        state.serialize_field("G->T", &self.mismatches().G.T)?;
        state.serialize_field("G->G", &self.mismatches().G.G)?;
        state.serialize_field("G->C", &self.mismatches().G.C)?;
        state.serialize_field("#C", &self.prednuc().C)?;
        state.serialize_field("C->A", &self.mismatches().C.A)?;
        state.serialize_field("C->T", &self.mismatches().C.T)?;
        state.serialize_field("C->G", &self.mismatches().C.G)?;
        state.serialize_field("C->C", &self.mismatches().C.C)?;
        state.end()
    }
}

// #[cfg(test)]
// mod test {
//     use bio_types::genome::Interval;
//     use bio_types::strand::Strand;
//
//     use crate::core::dna::NucCounts;
//     use crate::core::dna::Nucleotide;
//     use crate::core::mismatches::roi::{MismatchesSummary, MockROIMismatches};
//     use crate::core::workload::ROI;
//
//     use super::*;
//
//     #[test]
//     fn regions() {
//         let factory = |roi, strand, coverage, masked, sequence, ncounts| {
//             let mut dummy = MockROIMismatches::new();
//             let mut mismatches = MismatchesSummary::zeros();
//             mismatches.increment(sequence, ncounts);
//             let mut ncounts = NucCounts::zeros();
//             ncounts.increment(sequence);
//
//             dummy.expect_roi().return_const(roi);
//             dummy.expect_strand().return_const(strand);
//             dummy.expect_coverage().return_const(coverage);
//             dummy.expect_masked().return_const(masked);
//             dummy.expect_sequence().return_const(ncounts);
//             dummy.expect_mismatches().return_const(mismatches);
//             dummy
//         };
//
//         let workload = vec![
//             factory(
//                 ROI::new(Interval::new("1".into(), 1..20), "First".into(), vec![1..20]),
//                 Strand::Forward,
//                 12u32,
//                 0u32,
//                 &vec![Nucleotide::A, Nucleotide::Unknown],
//                 &vec![NucCounts::new(1, 10, 100, 1000), NucCounts::new(1, 1, 1, 1)],
//             ),
//             factory(
//                 ROI::new(Interval::new("3".into(), 30..32), "".into(), vec![30..31]),
//                 Strand::Unknown,
//                 190u32,
//                 1u32,
//                 &vec![Nucleotide::C, Nucleotide::G],
//                 &vec![NucCounts::new(1, 2, 3, 4), NucCounts::new(5, 6, 7, 8)],
//             ),
//         ];
//         let mut saveto = Vec::new();
//         super::regions(&mut saveto, workload);
//
//         let result = String::from_utf8(saveto).unwrap();
//         let expected = "chr\tstart\tend\tstrand\tname\tcoverage\t#masked\t\
//                               #A\t#T\t#G\t#C\t\
//                               A->A\tA->C\tA->G\tA->T\t\
//                               C->A\tC->C\tC->G\tC->T\t\
//                               G->A\tG->C\tG->G\tG->T\t\
//                               T->A\tT->C\tT->G\tT->T\n\
//                               1\t1\t20\t+\tFirst\t12\t0\t\
//                               1\t0\t0\t0\t\
//                               1\t10\t100\t1000\t\
//                               0\t0\t0\t0\t\
//                               0\t0\t0\t0\t\
//                               0\t0\t0\t0\n\
//                               3\t30\t32\t.\t\t190\t1\t\
//                               0\t0\t1\t1\t\
//                               0\t0\t0\t0\t\
//                               1\t2\t3\t4\t\
//                               5\t6\t7\t8\t\
//                               0\t0\t0\t0\n";
//         assert_eq!(&result, expected)
//     }
// }
