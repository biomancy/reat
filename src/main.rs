use std::sync::Arc;

use clap::{crate_authors, crate_name, crate_version, App, AppSettings, ArgMatches};
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use rayon::{ThreadPool, ThreadPoolBuilder};

use rada::cli;

const CREATE_THREAD_POOL_ERROR: &str = "Failed to initialize thread pool";

fn main() {
    let app = App::new(crate_name!())
        .author(crate_authors!("\n"))
        .version(crate_version!())
        .max_term_width(120)
        .setting(AppSettings::DeriveDisplayOrder)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("rois")
                .long_about("Quantify editing for the specified set of Regions Of Interest (ROIs)")
                .args(cli::roi::args()),
        )
        .subcommand(
            App::new("loci").long_about("Estimate editing per-loci for the whole genome.").args(cli::loci::args()),
        )
        .get_matches();

    let style = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner} {msg}")
        .tick_strings(&["▹▹▹▹▹", "▸▹▹▹▹", "▹▸▹▹▹", "▹▹▸▹▹", "▹▹▹▸▹", "▹▹▹▹▸", "▪▪▪▪▪"])
        .on_finish(ProgressFinish::AndLeave);
    let mbar = MultiProgress::new();
    let factory = || mbar.add(ProgressBar::new_spinner().with_style(style.clone()));

    let pbar = factory();
    pbar.set_message("Running...");

    ThreadPoolBuilder::new().num_threads(2).build().expect(CREATE_THREAD_POOL_ERROR).scope(|s| {
        // Render progress bar in the first thread
        s.spawn(|_| mbar.join().expect("Failed to render progress bar"));

        // Run subcommand in the second command
        s.spawn(|_| {
            match app.subcommand() {
                Some(("rois", matches)) => {
                    let core = cli::shared::args::CoreArgs::new(matches, factory);
                    ThreadPoolBuilder::new().num_threads(core.threads).build_global().expect(CREATE_THREAD_POOL_ERROR);
                    cli::roi::run(matches, core, factory);
                }
                Some(("loci", matches)) => {
                    let core = cli::shared::args::CoreArgs::new(matches, factory);
                    ThreadPoolBuilder::new().num_threads(core.threads).build_global().expect(CREATE_THREAD_POOL_ERROR);
                    cli::loci::run(matches, core, factory);
                }
                _ => panic!("Subcommand is not specified."),
            };

            // After finishing,
            pbar.finish_with_message("Finished!");
        });
    });
}
