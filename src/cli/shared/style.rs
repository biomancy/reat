use indicatif::{ProgressFinish, ProgressStyle};

pub fn running() -> ProgressStyle {
    ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-")
        .on_finish(ProgressFinish::AndLeave)
}

pub fn finished() -> ProgressStyle {
    ProgressStyle::default_bar().template("[{elapsed_precise}] {msg}").on_finish(ProgressFinish::AndLeave)
}
