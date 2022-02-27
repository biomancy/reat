use bio_types::strand::Strand;
use derive_getters::Dissolve;

use crate::core::mismatches::BatchedMismatches;

#[derive(Dissolve)]
pub struct StrandingContext<T: BatchedMismatches> {
    unstranded: Vec<T>,
    stranded: Vec<T>,
}

pub enum StrandingAlgoResult {
    AllElements(Strand),
    EachElement(Vec<Strand>),
}

impl<T: BatchedMismatches> StrandingContext<T> {
    pub fn new(items: Vec<T>) -> Self {
        let (unstranded, stranded) = items.into_iter().partition(|x| x.strand().is_unknown());
        Self { unstranded, stranded }
    }

    pub fn apply(&mut self, algo: impl Fn(&T) -> StrandingAlgoResult) {
        self.unstranded = std::mem::take(&mut self.unstranded)
            .into_iter()
            .filter_map(|mut x| match algo(&x) {
                StrandingAlgoResult::AllElements(strand) => match strand {
                    Strand::Forward | Strand::Reverse => {
                        x.restrand_all(strand);
                        self.stranded.push(x);
                        None
                    }
                    Strand::Unknown => Some(x),
                },
                StrandingAlgoResult::EachElement(strands) => {
                    let x = x.restrand(strands);
                    for elem in [x.forward, x.reverse] {
                        if let Some(e) = elem {
                            self.stranded.push(e);
                        }
                    }
                    x.unknown
                }
            })
            .collect();
    }
}
