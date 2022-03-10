use bio_types::strand::Strand;

use crate::core::mismatches::StrandingCounts;
use crate::core::strandutil::StrandedData;

pub fn maskvec<T>(mut v: Vec<T>, mask: &[bool]) -> Vec<T> {
    debug_assert_eq!(v.len(), mask.len());

    let mut index = 0;
    v.retain(|_| {
        index += 1;
        mask[index - 1]
    });
    v
}

pub fn select_strands<T>(v: Vec<T>, strands: &[Strand], cnts: &StrandingCounts) -> StrandedData<Vec<T>> {
    debug_assert_eq!(strands.len(), v.len());
    debug_assert_eq!(strands.len(), cnts.forward + cnts.reverse + cnts.unknown);

    let mut items = StrandedData {
        forward: Vec::with_capacity(cnts.forward),
        reverse: Vec::with_capacity(cnts.reverse),
        unknown: Vec::with_capacity(cnts.unknown),
    };

    for (v, &s) in v.into_iter().zip(strands) {
        items[s].push(v);
    }
    debug_assert_eq!(cnts.forward, items.forward.len());
    debug_assert_eq!(cnts.reverse, items.reverse.len());
    debug_assert_eq!(cnts.unknown, items.unknown.len());
    items
}

// TODO: units for these
