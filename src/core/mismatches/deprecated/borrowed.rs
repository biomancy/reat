impl SeparableIntervalMismatches for RefIntervalMismatches<'_> {
    fn split(self, indices: &[usize]) -> Vec<Self> {
        if indices.is_empty() {
            return vec![self];
        }

        debug_assert!(indices[0] > 0);
        let mut result = Vec::with_capacity(indices.len());

        let mut curind = 0;
        let (contig, start) = (self.interval.contig(), self.interval.range().start as u64);
        for &nextind in indices {
            result.push(Self {
                interval: Interval::new(contig.to_owned(), start + curind as u64..start + nextind as u64),
                strand: self.strand,
                refnuc: &self.refnuc[curind..nextind],
                ncounts: &self.ncounts[curind..nextind],
            });
            curind = nextind;
        }
        if *indices.last().unwrap() != self.refnuc.len() {
            let nextind = self.refnuc.len();
            result.push(Self {
                interval: Interval::new(contig.to_owned(), start + curind as u64..start + nextind as u64),
                strand: self.strand,
                refnuc: &self.refnuc[curind..],
                ncounts: &self.ncounts[curind..],
            });
        }

        debug_assert_eq!(result.len(), indices.len() + 1);
        result
    }

    fn extract(self, ranges: &[Range<usize>], into: &mut Vec<Self>) {
        let contig = self.interval.contig();
        let start = self.interval.range().start;

        into.extend(ranges.into_iter().map(|x| Self {
            interval: Interval::new(contig.to_owned(), start + x.start as u64..start + x.end as u64),
            strand: self.strand,
            refnuc: &self.refnuc[x.to_owned()],
            ncounts: &self.ncounts[x.to_owned()],
        }))
    }
}
