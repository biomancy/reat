use std::io::Write;

use bio_types::genome::{AbstractInterval, AbstractLocus};

use crate::core::stats::ROIBasedStat;
use crate::core::summary::{LocusSummary, ROISummary};

const OUTPUT_IO_ERROR: &str = "Failed to write loci summary to the output TSV file.";

pub fn loci(saveto: &mut impl Write, summary: Vec<LocusSummary>) {
    writeln!(saveto, "chr\tposition\tstrand\treference\tA\tC\tG\tT").expect(OUTPUT_IO_ERROR);

    for e in summary {
        writeln!(
            saveto,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.locus.contig(),
            e.locus.pos(),
            e.strand.strand_symbol(),
            e.refnuc,
            e.sequenced.A,
            e.sequenced.C,
            e.sequenced.G,
            e.sequenced.T
        )
        .expect(OUTPUT_IO_ERROR);
    }
}

#[cfg(test)]
mod test {
    use bio_types::genome::{Interval, Locus};
    use bio_types::strand::Strand;

    use crate::core::counting::NucCounts;
    use crate::core::dna::Nucleotide;

    use super::*;

    #[test]
    fn loci() {
        let dummy = vec![
            LocusSummary::new(Locus::new("1".into(), 10), Strand::Unknown, Nucleotide::Unknown, NucCounts::zeros()),
            LocusSummary::new(Locus::new("2".into(), 20), Strand::Forward, Nucleotide::G, NucCounts::new(1, 2, 3, 4)),
        ];
        let mut saveto = Vec::new();
        super::loci(&mut saveto, dummy);

        let result = String::from_utf8(saveto).unwrap();
        let expected = "chr\tposition\tstrand\treference\tA\tC\tG\tT\n\
                              1\t10\t.\tN\t0\t0\t0\t0\n\
                              2\t20\t+\tG\t1\t2\t3\t4\n";
        assert_eq!(&result, expected)
    }
}
