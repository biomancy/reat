use std::any::Any;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io;

use indicatif::ProgressBar;
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use serde::ser::Serialize;

use crate::cli::shared;
use crate::cli::shared::thread_cache::ThreadCache;
use crate::core::hooks::stats::EditingStatType;
use crate::core::hooks::stats::ROIEditingIndex;
use crate::core::mismatches::MismatchesVec;
use crate::core::runner::Runner;

const OUTPUT_IO_ERROR: &str = "Failed to write results to the output TSV file.";
const STATS_IO_ERROR: &str = "Failed to write statistics to the output TSV file.";

pub fn run<R, Edits, BM: MismatchesVec, Workload, W: io::Write>(
    workload: Vec<Workload>,
    runner: R,
    pbar: ProgressBar,
    saveto: &mut csv::Writer<W>,
    mut statsto: HashMap<EditingStatType, csv::Writer<W>>,
) -> csv::Result<()>
where
    Edits: Send + Serialize,
    Workload: Sized + Send,
    R: for<'runner> Runner<'runner, BM, Result = Vec<Edits>, Workload = Workload> + Clone + Send,
{
    // Callbacks to track progress
    pbar.set_style(shared::style::run::running());
    pbar.set_draw_delta((rayon::current_num_threads() * 10) as u64);
    pbar.set_length(workload.len() as u64);

    let ctxstore = ThreadCache::new(move || RefCell::new(runner.clone()));
    let edits: Vec<Edits> = workload
        .into_par_iter()
        .map(|w| {
            let result = ctxstore.get().borrow_mut().run(w);
            pbar.inc(1);
            result
        })
        .flatten()
        .collect();

    // Report the result
    pbar.set_style(shared::style::run::finished());
    pbar.finish_with_message(format!("Finished with {} items, total processed reads: TODO", edits.len()));

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
        saveto.serialize(edit).expect(STATS_IO_ERROR);
    }
    Ok(Default::default())
}
