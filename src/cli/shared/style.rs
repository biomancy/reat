use indicatif::{ProgressFinish, ProgressStyle};

pub mod parse {
    use super::*;

    pub fn with_progress() -> ProgressStyle {
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {spinner} {msg}")
            .tick_strings(&["▹▹▹▹▹", "▸▹▹▹▹", "▹▸▹▹▹", "▹▹▸▹▹", "▹▹▹▸▹", "▹▹▹▹▸", "▪▪▪▪▪"])
            .on_finish(ProgressFinish::AndLeave)
    }
}

pub mod run {
    use super::*;

    pub fn running() -> ProgressStyle {
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:60.cyan/blue} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-")
            .on_finish(ProgressFinish::AndLeave)
    }

    pub fn finished() -> ProgressStyle {
        ProgressStyle::default_bar().template("[{elapsed_precise}] {msg}").on_finish(ProgressFinish::AndLeave)
    }
}
