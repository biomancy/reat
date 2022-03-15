use std::ops::Range;

use bio_types::genome::{AbstractInterval, Interval, Position};
use derive_getters::{Dissolve, Getters};

use super::utils;

#[derive(Clone, PartialEq, Debug, Getters, Dissolve)]
pub struct SiteWorkload {
    interval: Interval,
    include: Vec<Range<u64>>,
}

impl AbstractInterval for SiteWorkload {
    fn contig(&self) -> &str {
        self.interval.contig()
    }

    fn range(&self) -> Range<Position> {
        self.interval.range()
    }
}

impl SiteWorkload {
    pub fn from_intervals(
        mut intervals: Vec<Interval>,
        binsize: u64,
        exclude: Option<Vec<impl AbstractInterval>>,
    ) -> Vec<SiteWorkload> {
        assert!(binsize > 0, "Binsize must be > 0");
        // Subtract excluded if needed
        if let Some(excluded) = exclude {
            intervals = utils::subtract(intervals, excluded)
                .into_iter()
                .flat_map(|x| {
                    let contig = x.inner.contig().to_owned();
                    x.retained.into_iter().map(move |piece| Interval::new(contig.clone(), piece))
                })
                .collect();
        }

        // Bin and transform to the workload
        utils::split(intervals, binsize)
            .into_iter()
            .map(|x| {
                let include = vec![x.range().start..x.range().end];
                SiteWorkload { interval: x, include }
            })
            .collect()
    }
}
