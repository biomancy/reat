use std::io::Write;

use bio_types::genome::{AbstractInterval, AbstractLocus};

use crate::rada::summary::{IntervalSummary, LocusSummary};

pub fn loci(mut saveto: impl Write, summary: Vec<LocusSummary>) {
    write!(
        saveto,
        "chr\tposition\tstrand\t\
        reference\t#A\t#C\t#G\t#T\n"
    );

    for e in summary {
        write!(saveto, "{}\t{}\t{}\t", e.locus.contig(), e.locus.pos(), e.strand,);
        write!(saveto, "{}\t{}\t{}\t{}\t{}\n", e.refnuc, e.sequenced.A, e.sequenced.C, e.sequenced.G, e.sequenced.T);
    }
}

pub fn regions(mut saveto: impl Write, summary: Vec<IntervalSummary>) {
    write!(
        saveto,
        "chr\tstart\tend\tstrand\tname\t\
        A->A\tA->C\tA->G\tA->T\t\"
        C->A\tC->C\tC->G\tC->T\t\"
        G->A\tG->C\tG->G\tG->T\t\"
        T->A\tT->C\tT->G\tT->T\t\n"
    );

    for e in summary {
        let (contig, range) = (e.interval.contig(), e.interval.range());
        write!(saveto, "{}\t{}\t{}\t{}\t{}\t", contig, range.start, range.end, e.strand.strand_symbol(), e.name);
        let m = e.mismatches;
        write!(saveto, "{}\t{}\t{}\t{}\t", m.A.A, m.A.C, m.A.G, m.A.T);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.C.A, m.C.C, m.C.G, m.C.T);
        write!(saveto, "{}\t{}\t{}\t{}\t", m.G.A, m.G.C, m.G.G, m.G.T);
        write!(saveto, "{}\t{}\t{}\t{}\n", m.T.A, m.T.C, m.T.G, m.T.T);
    }
}
