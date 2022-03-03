use std::cmp::Ordering;
use std::cmp::{max, min};
use std::ops::Range;

use bio_types::genome::AbstractInterval;

#[derive(Eq, PartialEq, Debug)]
pub struct MaskedInterval<T: AbstractInterval> {
    pub inner: T,
    pub retained: Vec<Range<u64>>,
}

fn cmp(this: &impl AbstractInterval, other: &impl AbstractInterval) -> Ordering {
    let from_contig = this.contig().cmp(other.contig());
    if !from_contig.is_eq() {
        return from_contig;
    } else {
        this.range().start.cmp(&other.range().start)
    }
}

fn chop_pieces(pieces: &mut Vec<Range<u64>>, by: &Range<u64>, buffer: &mut Vec<Range<u64>>) {
    buffer.clear();

    for piece in pieces.iter() {
        // Outside of the range
        if piece.start >= by.end || piece.end <= by.start {
            buffer.push(piece.start..piece.end);
            continue;
        }
        // Complete overlap
        if piece.start >= by.start && piece.end <= by.end {
            continue;
        }
        // Partial overlap
        let overlap = max(piece.start, by.start)..min(piece.end, by.end);
        debug_assert!(overlap.end - overlap.start > 0);
        if overlap.start > piece.start {
            buffer.push(piece.start..overlap.start)
        }
        if overlap.end < piece.end {
            buffer.push(overlap.end..piece.end)
        }
    }
    std::mem::swap(pieces, buffer);
}

fn process_window<T: AbstractInterval>(
    mut window: Vec<MaskedInterval<T>>,
    subtract: &impl AbstractInterval,
    saveto: &mut Vec<MaskedInterval<T>>,
    windowbuf: &mut Vec<MaskedInterval<T>>,
    piecesbuf: &mut Vec<Range<u64>>,
) -> Vec<MaskedInterval<T>> {
    windowbuf.clear();

    for mut w in window.drain(..).into_iter() {
        // Behind the subtract => save to the result
        if w.inner.contig() < subtract.contig() || w.inner.range().end <= subtract.range().start {
            if !w.retained.is_empty() {
                saveto.push(w);
            }
            continue;
        }
        // Ahead of the subtract => save to the buffer
        if w.inner.contig() > subtract.contig() || w.inner.range().start >= subtract.range().end {
            windowbuf.push(w);
            continue;
        }
        // Must overlap
        debug_assert!(w.inner.contig() == subtract.contig());
        debug_assert!(
            min(w.inner.range().end, subtract.range().end) - max(w.inner.range().start, subtract.range().start) > 0
        );
        chop_pieces(&mut w.retained, &subtract.range(), piecesbuf);
        if !w.retained.is_empty() {
            windowbuf.push(w);
        }
    }
    std::mem::swap(windowbuf, &mut window);
    window
}

pub fn subtract<T: AbstractInterval>(
    mut inters: Vec<T>,
    mut subtract: Vec<impl AbstractInterval>,
) -> Vec<MaskedInterval<T>> {
    // Sort intervals by coordinate
    inters.sort_by(cmp);
    subtract.sort_by(cmp);

    let mut saveto = Vec::with_capacity(inters.len());
    let mut window = Vec::with_capacity(100);

    let mut piecesbuf = Vec::with_capacity(100);
    let mut windowbuf = Vec::with_capacity(100);

    let (mut interit, mut subit) = (inters.into_iter(), subtract.into_iter());
    let (mut nextiter, mut nextsub) = (interit.next(), subit.next());

    loop {
        match (nextiter, nextsub) {
            (Some(i), Some(s)) => {
                // Skip intervals before the subtract
                if i.contig() < s.contig() && i.range().end <= s.range().start {
                    let range = i.range();
                    saveto.push(MaskedInterval { inner: i, retained: vec![range] });
                    nextiter = interit.next();
                    nextsub = Some(s);
                    continue;
                }

                // Fill the window while possible
                if i.contig() == s.contig() && (i.range().start < s.range().end) {
                    let range = i.range();
                    window.push(MaskedInterval { inner: i, retained: vec![range] });
                    nextiter = interit.next();
                    nextsub = Some(s);
                    continue;
                }

                // Subtract if possible
                window = process_window(window, &s, &mut saveto, &mut windowbuf, &mut piecesbuf);

                nextiter = Some(i);
                nextsub = subit.next();
            }
            (Some(i), None) => {
                // Nothing to subtract, simply extend the result
                let range = i.range();
                saveto.push(MaskedInterval { inner: i, retained: vec![range] });
                // Update values
                nextiter = interit.next();
                nextsub = None;
            }
            (None, Some(s)) => {
                window = process_window(window, &s, &mut saveto, &mut windowbuf, &mut piecesbuf);
                // Update values
                nextiter = None;
                nextsub = subit.next();
            }
            (None, None) => {
                window.retain(|x| !x.retained.is_empty());
                if !window.is_empty() {
                    saveto.extend(window.into_iter());
                }
                break;
            }
        }
    }
    saveto
}

#[cfg(test)]
mod tests {
    use bio_types::genome::Interval;
    use itertools::Itertools;

    use super::*;

    fn mwork(chrom: &str, ranges: Vec<Range<u64>>) -> Vec<Interval> {
        ranges.into_iter().map(|x| Interval::new(chrom.into(), x)).collect()
    }

    #[test]
    fn simple() {
        let inter = mwork("1", vec![2..8]);
        for (sub, expected) in
            [(4..6, vec![2..4, 6..8]), (0..5, vec![5..8]), (6..9, vec![2..6]), (0..2, vec![2..8]), (8..10, vec![2..8])]
        {
            let sub = mwork("1", vec![sub]);
            let expected = vec![MaskedInterval { inner: inter[0].clone(), retained: expected }];
            let result = subtract(inter.clone(), sub);
            assert_eq!(expected, result);
        }
    }

    #[test]
    fn empty() {
        let inter = mwork("1", vec![2..8]);
        // Empty subtract
        let expected = vec![MaskedInterval { inner: inter[0].clone(), retained: vec![2..8] }];
        let result = subtract(inter.clone(), Vec::<Interval>::new());
        assert_eq!(expected, result);

        // Empty inter
        let sub = mwork("1", vec![0..100]);
        let result = subtract(Vec::<Interval>::new(), sub);
        assert!(result.is_empty());

        // Empty result
        for sub in [vec![2..8], vec![0..10], vec![0..2, 2..4, 4..7, 7..8]] {
            let sub = mwork("1", sub);
            let result = subtract(inter.clone(), sub);
            assert!(result.is_empty());
        }
    }

    #[test]
    fn complex() {
        // Test case 1
        let inter_1 = mwork("1", vec![1..3, 2..8, 2..4, 5..6, 7..9, 10..11, 0..13]);
        let sub_1 = mwork("1", vec![2..3, 2..5, 5..9, 7..12]);
        let expect_1 = vec![
            MaskedInterval { inner: inter_1[0].clone(), retained: vec![1..2] },
            MaskedInterval { inner: inter_1[6].clone(), retained: vec![0..2, 12..13] },
        ];
        assert_eq!(subtract(inter_1.clone(), sub_1.clone()), expect_1);
        // Test case 2
        let inter_2 = mwork("2", vec![0..2, 4..12, 0..3, 4..5, 6..9, 0..1, 4..5, 6..8, 9..13, 6..7]);
        let sub_2 = mwork("2", vec![0..1, 0..1, 2..3, 3..4, 5..6, 7..10, 11..13, 1..4, 5..7, 12..14]);
        let expect_2 = vec![
            MaskedInterval { inner: inter_2[3].clone(), retained: vec![4..5] },
            MaskedInterval { inner: inter_2[6].clone(), retained: vec![4..5] },
            MaskedInterval { inner: inter_2[1].clone(), retained: vec![4..5, 10..11] },
            MaskedInterval { inner: inter_2[8].clone(), retained: vec![10..11] },
        ];
        assert_eq!(subtract(inter_2.clone(), sub_2.clone()), expect_2);
        // Test case 3
        let inter_3 = mwork("3", vec![0..2, 1..3, 2..3, 6..7, 6..8, 6..10, 8..11]);
        let sub_3 = mwork("3", vec![1..2, 4..5, 4..6, 4..6, 7..10, 12..15]);
        let expect_3 = vec![
            MaskedInterval { inner: inter_3[0].clone(), retained: vec![0..1] },
            MaskedInterval { inner: inter_3[1].clone(), retained: vec![2..3] },
            MaskedInterval { inner: inter_3[2].clone(), retained: vec![2..3] },
            MaskedInterval { inner: inter_3[3].clone(), retained: vec![6..7] },
            MaskedInterval { inner: inter_3[4].clone(), retained: vec![6..7] },
            MaskedInterval { inner: inter_3[5].clone(), retained: vec![6..7] },
            MaskedInterval { inner: inter_3[6].clone(), retained: vec![10..11] },
        ];
        assert_eq!(subtract(inter_3.clone(), sub_3.clone()), expect_3);
        // Test case = sum of the above
        let inter = [inter_1, inter_2, inter_3].into_iter().flatten().collect_vec();
        let sub = [sub_1, sub_2, sub_3].into_iter().flatten().collect_vec();
        let expect = [expect_1, expect_2, expect_3].into_iter().flatten().collect_vec();

        assert_eq!(subtract(inter.clone(), sub.clone()), expect);
    }
}
