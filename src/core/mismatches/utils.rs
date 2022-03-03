use bio_types::strand::Strand;

use crate::core::mismatches::StrandingCounts;

pub fn maskvec<T>(mut v: Vec<T>, mask: &[bool]) -> Vec<T> {
    debug_assert_eq!(v.len(), mask.len());

    let mut index = 0;
    v.retain(|_| {
        index += 1;
        mask[index - 1]
    });
    v
}

pub fn select_strands<T>(v: Vec<T>, strands: &[Strand], cnts: &StrandingCounts) -> (Vec<T>, Vec<T>, Vec<T>) {
    debug_assert_eq!(strands.len(), v.len());
    debug_assert_eq!(strands.len(), cnts.forward + cnts.reverse + cnts.unknown);

    let (mut fwd, mut rev, mut unk) =
        (Vec::with_capacity(cnts.forward), Vec::with_capacity(cnts.reverse), Vec::with_capacity(cnts.unknown));

    for (v, s) in v.into_iter().zip(strands) {
        match s {
            Strand::Forward => fwd.push(v),
            Strand::Reverse => rev.push(v),
            Strand::Unknown => unk.push(v),
        }
    }
    (fwd, rev, unk)
}
