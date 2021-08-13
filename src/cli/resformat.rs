use std::io::Write;

use bio_types::genome::{AbstractInterval, AbstractLocus};

use crate::rada::summary::{IntervalSummary, LocusSummary};

const IO_ERROR: &str = "Failed to write to the output file.";

pub fn loci(mut saveto: impl Write, summary: Vec<LocusSummary>) {
    writeln!(saveto, "chr\tposition\tstrand\treference\tA\tC\tG\tT").expect(IO_ERROR);

    for e in summary {
        writeln!(
            saveto,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.locus.contig(),
            e.locus.pos(),
            e.strand,
            e.refnuc,
            e.sequenced.A,
            e.sequenced.C,
            e.sequenced.G,
            e.sequenced.T
        )
        .expect(IO_ERROR);
    }
}

pub fn regions(mut saveto: impl Write, summary: Vec<IntervalSummary>) {
    writeln!(saveto, "chr\tstart\tend\tstrand\tname\tA->A\tA->C\tA->G\tA->T\tC->A\tC->C\tC->G\tC->T\tG->A\tG->C\tG->G\tG->T\tT->A\tT->C\tT->G\tT->T\t")
    .expect(IO_ERROR);

    for e in summary {
        let (contig, range, strand) = (e.interval.contig(), e.interval.range(), e.strand.strand_symbol());
        write!(saveto, "{}\t{}\t{}\t{}\t{}\t", contig, range.start, range.end, strand, e.name).expect(IO_ERROR);
        let m = e.mismatches;
        write!(saveto, "{}\t{}\t{}\t{}\t", m.A.A, m.A.C, m.A.G, m.A.T).expect(IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.C.A, m.C.C, m.C.G, m.C.T).expect(IO_ERROR);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.G.A, m.G.C, m.G.G, m.G.T).expect(IO_ERROR);
        writeln!(saveto, "{}\t{}\t{}\t{}", m.T.A, m.T.C, m.T.G, m.T.T).expect(IO_ERROR);
    }
}
