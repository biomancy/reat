use indicatif::{ProgressFinish, ProgressStyle};

pub fn pbar() -> ProgressStyle {
    ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-")
        .on_finish(ProgressFinish::AndLeave)
}
