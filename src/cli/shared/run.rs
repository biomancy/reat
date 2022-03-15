use std::any::Any;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io;

use bio_types::strand::Strand;
use indicatif::ProgressBar;
use itertools::min;
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use serde::ser::Serialize;

use crate::cli::shared;
use crate::cli::shared::thread_cache::ThreadCache;
use crate::core::hooks::stats::EditingStatType;
use crate::core::hooks::stats::ROIEditingIndex;
use crate::core::mismatches::{Batch, MismatchesVec};
use crate::core::runner::Runner;
use crate::core::strandutil::Stranded;

const OUTPUT_IO_ERROR: &str = "Failed to write results to the output TSV file.";
const STATS_IO_ERROR: &str = "Failed to write statistics to the output TSV file.";

pub fn run<RunnerT, Mismatches, Workload, W: io::Write>(
    workload: Vec<Workload>,
    runner: RunnerT,
    pbar: ProgressBar,
    saveto: &mut csv::Writer<W>,
    mut statsto: HashMap<EditingStatType, csv::Writer<W>>,
) -> csv::Result<()>
where
    Mismatches: Send + MismatchesVec,
    Workload: Sized + Send,
    RunnerT: for<'runner> Runner<'runner, Mismatches, Workload = Workload> + Clone + Send,
{
    // Callbacks to track progress
    pbar.set_style(shared::style::run::running());

    // let delta = std::cmp::min(rayon::current_num_threads() * 10, workload.len() / 10 + 1);
    // pbar.set_draw_delta(delta as u64);

    pbar.set_length(workload.len() as u64);

    let ctxstore = ThreadCache::new(move || RefCell::new(runner.clone()));
    let edits: Vec<Batch<Mismatches>> = workload
        // .into_par_iter()
        .into_iter()
        .filter_map(|w| {
            let mut result = ctxstore.get().borrow_mut().run(w);
            pbar.inc(1);
            result
        })
        .collect();

    // Report the result
    pbar.set_style(shared::style::run::finished());
    let reads: Stranded<u32> = edits.iter().map(|x| x.mapped).fold(Default::default(), |a, b| a + b);
    pbar.finish_with_message(format!("Finished with {} items, processed reads: {}", edits.len(), reads));

    // Group stats by type
    let stats = ctxstore.dissolve().map(|x| x.into_inner().stats()).flatten();
    let mut grouped: HashMap<EditingStatType, Vec<Box<dyn Any>>> = HashMap::new();
    for stat in stats {
        let (typed, any) = stat.into_any();
        grouped.entry(typed).or_default().push(any);
    }

    // Collapse identical stats & write them into requested serializers
    for (k, v) in grouped.into_iter() {
        if let Some(serializer) = statsto.get_mut(&k) {
            match k {
                EditingStatType::ROIEditingIndex => {
                    serializer.serialize(ROIEditingIndex::collapse(v)).expect(OUTPUT_IO_ERROR)
                }
            };
        };
    }

    // Serialize edits
    for edit in edits {
        for x in [edit.items, edit.retained] {
            x.forward.to_csv(saveto).expect(OUTPUT_IO_ERROR);
            x.reverse.to_csv(saveto).expect(OUTPUT_IO_ERROR);
            x.unknown.to_csv(saveto).expect(OUTPUT_IO_ERROR);
        }
    }
    Ok(())
}
