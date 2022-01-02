use std::io::Write;

use crate::core::mismatches::site::SiteMismatches;
use bio_types::genome::AbstractLocus;

const OUTPUT_IO_ERROR: &str = "Failed to write site filters to the output TSV file.";

pub fn sites(saveto: &mut impl Write, summary: Vec<impl SiteMismatches>) {
    writeln!(saveto, "chr\tposition\tstrand\treference\tA\tC\tG\tT").expect(OUTPUT_IO_ERROR);

    for e in summary {
        writeln!(
            saveto,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.locus().contig(),
            e.locus().pos(),
            e.strand().strand_symbol(),
            e.refnuc(),
            e.ncounts().A,
            e.ncounts().C,
            e.ncounts().G,
            e.ncounts().T
        )
        .expect(OUTPUT_IO_ERROR);
    }
}

#[cfg(test)]
mod test {
    use bio_types::genome::Locus;
    use bio_types::strand::Strand;

    use crate::core::dna::NucCounts;
    use crate::core::dna::Nucleotide;

    use super::*;
    use crate::core::mismatches::site::MockSiteMismatches;

    #[test]
    fn loci() {
        let factory = |locus, strand, nuc, counts| {
            let mut dummy = MockSiteMismatches::new();
            dummy.expect_locus().return_const(locus);
            dummy.expect_strand().return_const(strand);
            dummy.expect_refnuc().return_const(nuc);
            dummy.expect_ncounts().return_const(counts);
            dummy
        };

        let workload = vec![
            factory(Locus::new("1".into(), 10), Strand::Unknown, Nucleotide::Unknown, NucCounts::zeros()),
            factory(Locus::new("2".into(), 20), Strand::Forward, Nucleotide::G, NucCounts::new(1, 2, 3, 4)),
        ];
        let mut saveto = Vec::new();
        super::sites(&mut saveto, workload);

        let result = String::from_utf8(saveto).unwrap();
        let expected = "chr\tposition\tstrand\treference\tA\tC\tG\tT\n\
                              1\t10\t.\tN\t0\t0\t0\t0\n\
                              2\t20\t+\tG\t1\t2\t3\t4\n";
        assert_eq!(&result, expected)
    }
}
