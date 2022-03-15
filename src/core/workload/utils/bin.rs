use std::cmp::Ordering;

use bio_types::genome::{AbstractInterval, Interval};

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Bin<T> {
    pub bin: Interval,
    pub items: Vec<T>,
}

fn cmp(this: &impl AbstractInterval, other: &impl AbstractInterval) -> Ordering {
    let from_contig = this.contig().cmp(other.contig());
    if !from_contig.is_eq() {
        return from_contig;
    }

    let from_start = this.range().start.cmp(&other.range().start);
    if !from_start.is_eq() {
        return from_start;
    }

    let len = this.range().end - this.range().start;
    let otherlen = other.range().end - other.range().start;

    // Intervals with larget size must go first!
    len.cmp(&otherlen).reverse()
}

fn infer_bin<T: AbstractInterval>(seed: &T, binsize: u64) -> Interval {
    let end = seed.range().end.max(seed.range().start + binsize);
    Interval::new(seed.contig().to_string(), seed.range().start..end)
}

pub fn bin<T: AbstractInterval>(mut workloads: Vec<T>, maxbinsize: u64) -> Vec<Bin<T>> {
    if workloads.is_empty() {
        return vec![];
    }
    workloads.sort_by(cmp);

    let mut result = Vec::with_capacity(workloads.len());

    let mut iter = workloads.into_iter();
    let first = iter.next().unwrap();
    let mut bin = infer_bin(&first, maxbinsize);
    let mut buffer = vec![first];
    let mut maxend = buffer[0].range().end;

    for work in iter {
        // outside of the bin
        if work.contig() != bin.contig() || work.range().end > bin.range().end {
            debug_assert!(!buffer.is_empty());
            if maxend < bin.range().end {
                bin.range_mut().end = maxend
            }

            // Save current
            result.push(Bin { bin, items: buffer });
            // Start the new one
            bin = infer_bin(&work, maxbinsize);
            buffer = vec![work];
            maxend = buffer[0].range().end;
        } else {
            maxend = maxend.max(work.range().end);
            buffer.push(work)
        }
    }

    // Save results for the last bin
    if maxend < bin.range().end {
        bin.range_mut().end = maxend
    }
    result.push(Bin { bin, items: buffer });

    result
}

fn split_interval(chromosome: &str, len: u64, binsize: u64) -> Vec<Interval> {
    let total_bins = (len + binsize - 1) / binsize;
    let mut bins: Vec<Interval> = Vec::with_capacity(total_bins as usize);

    for bin in 1..(total_bins + 1) {
        let start = (bin - 1) * binsize;
        let end = std::cmp::min(bin * binsize, len);
        debug_assert!(end - start <= binsize);

        bins.push(Interval::new(chromosome.to_owned(), start..end))
    }
    bins
}

pub fn split(intervals: Vec<impl AbstractInterval>, binsize: u64) -> Vec<Interval> {
    intervals.into_iter().flat_map(|x| split_interval(x.contig(), x.range().end - x.range().start, binsize)).collect()
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use bio_types::genome::Position;

    use super::*;

    fn interval(chr: &str, range: Range<Position>) -> Interval {
        Interval::new(chr.to_string(), range)
    }

    fn validate(input: Vec<Interval>, expected: Vec<Bin<Interval>>, binsizes: &[u64]) {
        for binsize in binsizes {
            let result = bin(input.clone(), *binsize);
            assert_eq!(result, expected);
        }
    }

    #[test]
    fn non_overlapping() {
        let inter = vec![
            interval("chr1", 10..20),
            interval("chr1", 30..40),
            interval("chr1", 50..60),
            interval("chr1", 70..80),
        ];

        let expected = vec![
            Bin { bin: interval("chr1", 10..20), items: vec![inter[0].clone()] },
            Bin { bin: interval("chr1", 30..40), items: vec![inter[1].clone()] },
            Bin { bin: interval("chr1", 50..60), items: vec![inter[2].clone()] },
            Bin { bin: interval("chr1", 70..80), items: vec![inter[3].clone()] },
        ];
        validate(inter.clone(), expected, &[1, 15, 29]);

        let expected = vec![
            Bin { bin: interval("chr1", 10..40), items: inter[0..2].to_vec() },
            Bin { bin: interval("chr1", 50..80), items: inter[2..].to_vec() },
        ];
        validate(inter.clone(), expected, &[30, 40, 49]);

        let expected = vec![Bin { bin: interval("chr1", 10..80), items: inter.clone() }];
        validate(inter, expected, &[70]);
    }

    #[test]
    fn overlapping() {
        let inter = vec![
            interval("1", 0..3),
            interval("1", 2..5),
            interval("1", 3..7),
            interval("1", 4..8),
            interval("1", 4..8),
        ];

        let expected = vec![
            Bin { bin: interval("1", 0..3), items: vec![inter[0].clone()] },
            Bin { bin: interval("1", 2..5), items: vec![inter[1].clone()] },
            Bin { bin: interval("1", 3..7), items: vec![inter[2].clone()] },
            Bin { bin: interval("1", 4..8), items: vec![inter[3].clone(), inter[4].clone()] },
        ];
        validate(inter.clone(), expected, &[1, 2, 4]);

        let expected = vec![
            Bin { bin: interval("1", 0..5), items: inter[0..2].to_vec() },
            Bin { bin: interval("1", 3..8), items: inter[2..].to_vec() },
        ];
        validate(inter.clone(), expected, &[5, 6]);

        let expected = vec![Bin { bin: interval("1", 0..8), items: inter.clone() }];
        validate(inter, expected, &[8, 100]);
    }

    #[test]
    fn none() {
        validate(vec![], vec![], &[1, 2, 3, 4]);
    }

    #[test]
    fn single() {
        let inter = interval("MT", 100..200);
        let expected = vec![Bin { bin: interval("MT", 100..200), items: vec![inter.clone()] }];

        validate(vec![inter], expected, &[1, 100, 500]);
    }

    #[test]
    fn complex() {
        let inter = vec![
            interval("1", 1..3),
            interval("1", 2..3),
            interval("1", 100..110),
            interval("2", 0..200),
            interval("2", 30..50),
            interval("2", 110..120),
            interval("3", 0..10),
            interval("3", 10..20),
            interval("3", 20..30),
            interval("3", 30..40),
            interval("4", 1..50),
            interval("4", 1..25),
            interval("4", 1..13),
            interval("4", 1..500),
        ];

        let chr4bin = [inter[13].clone(), inter[10].clone(), inter[11].clone(), inter[12].clone()].to_vec();
        let chr4bin = Bin { bin: interval("4", 1..500), items: chr4bin };

        let expected = vec![
            Bin { bin: interval("1", 1..110), items: inter[..3].to_vec() },
            Bin { bin: interval("2", 0..200), items: inter[3..6].to_vec() },
            Bin { bin: interval("3", 0..40), items: inter[6..10].to_vec() },
            chr4bin.clone(),
        ];
        validate(inter.clone(), expected, &[200, 300]);

        let expected = vec![
            Bin { bin: interval("1", 1..3), items: inter[..2].to_vec() },
            Bin { bin: interval("1", 100..110), items: inter[2..3].to_vec() },
            Bin { bin: interval("2", 0..200), items: inter[3..6].to_vec() },
            Bin { bin: interval("3", 0..40), items: inter[6..10].to_vec() },
            chr4bin,
        ];
        validate(inter, expected, &[50]);
    }

    #[test]
    fn split_interval() {
        let workload = interval("chr1", 0..284);
        let expected = vec![interval("chr1", 0..100), interval("chr1", 100..200), interval("chr1", 200..284)];

        assert_eq!(split(vec![workload], 100), expected);

        let result = split(vec![interval("2", 0..10)], 1000);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], interval("2", 0..10));
    }
}
