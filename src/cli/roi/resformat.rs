use std::io::Write;

use bio_types::genome::{AbstractInterval, AbstractLocus};

use crate::core::stats::ROIBasedStat;
use crate::core::summary::{LocusSummary, ROISummary};

const OUTPUT_IO_ERROR: &str = "Failed to write ROI summary to the output TSV file.";
const STATS_IO_ERROR: &str = "Failed to write statistics to the output TSV file.";

pub fn regions(saveto: &mut impl Write, summary: Vec<ROISummary>) {
    writeln!(saveto, "chr\tstart\tend\tstrand\tname\t#A\t#T\t#G\t#C\tA->A\tA->C\tA->G\tA->T\tC->A\tC->C\tC->G\tC->T\tG->A\tG->C\tG->G\tG->T\tT->A\tT->C\tT->G\tT->T")
        .expect(OUTPUT_IO_ERROR);

    for e in summary {
        let (contig, range, strand) = (e.interval.contig(), e.interval.range(), e.strand.strand_symbol());
        write!(saveto, "{}\t{}\t{}\t{}\t{}\t", contig, range.start, range.end, strand, e.name).expect(OUTPUT_IO_ERROR);
        let seq = e.sequenced;
        write!(saveto, "{}\t{}\t{}\t{}\t", seq.A, seq.T, seq.G, seq.C).expect(OUTPUT_IO_ERROR);
        let m = e.mismatches;
        write!(saveto, "{}\t{}\t{}\t{}\t", m.A.A, m.A.C, m.A.G, m.A.T).expect(OUTPUT_IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.C.A, m.C.C, m.C.G, m.C.T).expect(OUTPUT_IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.G.A, m.G.C, m.G.G, m.G.T).expect(OUTPUT_IO_ERROR);
        writeln!(saveto, "{}\t{}\t{}\t{}", m.T.A, m.T.C, m.T.G, m.T.T).expect(OUTPUT_IO_ERROR);
    }
}

pub fn statistic<T: ROIBasedStat>(files: &str, saveto: &mut impl Write, stat: &T, add_header: bool) {
    if add_header {
        writeln!(saveto, "files\t{}", T::header()).expect(STATS_IO_ERROR);
    }
    writeln!(saveto, "{}\t{}", files, stat.row()).expect(STATS_IO_ERROR);
}

#[cfg(test)]
mod test {
    use bio_types::genome::{Interval, Locus};
    use bio_types::strand::Strand;

    use crate::core::counting::NucCounts;
    use crate::core::dna::Nucleotide;

    use super::*;

    #[test]
    fn regions() {
        let dummy = vec![
            ROISummary::from_counts(
                Interval::new("1".into(), 1..20),
                "First".into(),
                Strand::Forward,
                &vec![Nucleotide::A, Nucleotide::Unknown],
                &vec![NucCounts::new(1, 10, 100, 1000), NucCounts::new(1, 1, 1, 1)],
            ),
            ROISummary::from_counts(
                Interval::new("3".into(), 30..31),
                "".into(),
                Strand::Unknown,
                &vec![Nucleotide::C, Nucleotide::G],
                &vec![NucCounts::new(1, 2, 3, 4), NucCounts::new(5, 6, 7, 8)],
            ),
        ];
        let mut saveto = Vec::new();
        super::regions(&mut saveto, dummy);

        let result = String::from_utf8(saveto).unwrap();
        let expected = "chr\tstart\tend\tstrand\tname\t\
                              #A\t#T\t#G\t#C\t\
                              A->A\tA->C\tA->G\tA->T\t\
                              C->A\tC->C\tC->G\tC->T\t\
                              G->A\tG->C\tG->G\tG->T\t\
                              T->A\tT->C\tT->G\tT->T\n\
                              1\t1\t20\t+\tFirst\t\
                              1\t0\t0\t0\t\
                              1\t10\t100\t1000\t\
                              0\t0\t0\t0\t\
                              0\t0\t0\t0\t\
                              0\t0\t0\t0\n\
                              3\t30\t31\t.\t\t\
                              0\t0\t1\t1\t\
                              0\t0\t0\t0\t\
                              1\t2\t3\t4\t\
                              5\t6\t7\t8\t\
                              0\t0\t0\t0\n";
        assert_eq!(&result, expected)
    }
}
