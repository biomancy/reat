use bio_types::strand::Strand;

use crate::core::mismatches::StrandingCounts;
use crate::core::stranding::predict::StrandingAlgoResult;

pub fn stranditer(strands: impl Iterator<Item = Strand>) -> StrandingAlgoResult {
    let mut counts = StrandingCounts::default();
    let strands = strands
        .map(|s| {
            counts[s] += 1;
            s
        })
        .collect();

    match (counts.forward > 0, counts.reverse > 0, counts.unknown > 0) {
        (true, false, false) => StrandingAlgoResult::AllElements(Strand::Forward),
        (false, true, false) => StrandingAlgoResult::AllElements(Strand::Reverse),
        (false, false, true) => StrandingAlgoResult::AllElements(Strand::Unknown),
        _ => StrandingAlgoResult::EachElement((strands, counts)),
    }
}
