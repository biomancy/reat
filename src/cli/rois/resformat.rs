use std::io::Write;

use bio_types::genome::AbstractInterval;

use crate::core::io::Table;
use crate::core::mismatches::roi::ROIMismatches;

const OUTPUT_IO_ERROR: &str = "Failed to write ROI filters to the output TSV file.";
const STATS_IO_ERROR: &str = "Failed to write statistics to the output TSV file.";

pub fn regions(saveto: &mut impl Write, summary: Vec<impl ROIMismatches>) {
    writeln!(saveto, "chr\tstart\tend\tstrand\tname\tcoverage\t#masked\t#A\t#T\t#G\t#C\tA->A\tA->C\tA->G\tA->T\tC->A\tC->C\tC->G\tC->T\tG->A\tG->C\tG->G\tG->T\tT->A\tT->C\tT->G\tT->T")
        .expect(OUTPUT_IO_ERROR);

    for e in summary {
        write!(
            saveto,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            e.roi().contig(),
            e.roi().original().range().start,
            e.roi().original().range().end,
            e.strand().strand_symbol(),
            e.roi().name(),
            e.coverage(),
            e.masked()
        )
        .expect(OUTPUT_IO_ERROR);
        let seq = e.prednuc();
        write!(saveto, "{}\t{}\t{}\t{}\t", seq.A, seq.T, seq.G, seq.C).expect(OUTPUT_IO_ERROR);
        let m = e.mismatches();
        write!(saveto, "{}\t{}\t{}\t{}\t", m.A.A, m.A.C, m.A.G, m.A.T).expect(OUTPUT_IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.C.A, m.C.C, m.C.G, m.C.T).expect(OUTPUT_IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.G.A, m.G.C, m.G.G, m.G.T).expect(OUTPUT_IO_ERROR);
        writeln!(saveto, "{}\t{}\t{}\t{}", m.T.A, m.T.C, m.T.G, m.T.T).expect(OUTPUT_IO_ERROR);
    }
}

pub fn statheader<T: Table, W: Write>(saveto: &mut W) {
    writeln!(saveto, "Run name\t{}", T::header().join("\t")).expect(STATS_IO_ERROR);
}

pub fn statistic<T: Table>(rname: &str, saveto: &mut impl Write, stat: &T) {
    writeln!(saveto, "{}\t{}", rname, stat.row().join("\t")).expect(STATS_IO_ERROR);
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
