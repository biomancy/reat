use clap::{crate_authors, crate_name, crate_version, App, AppSettings, ArgMatches};
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use rayon::ThreadPoolBuilder;

use rada::cli;
use rada::cli::shared::args::CoreArgs;

const CREATE_THREAD_POOL_ERROR: &str = "Failed to initialize thread pool";
const RENDER_PROGRESS_ERROR: &str = "Failed to render progress bar";

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

    // Setup progress tracking
    let style = ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {spinner} {msg}")
        .tick_strings(&["▹▹▹▹▹", "▸▹▹▹▹", "▹▸▹▹▹", "▹▹▸▹▹", "▹▹▹▸▹", "▹▹▹▹▸", "▪▪▪▪▪"])
        .on_finish(ProgressFinish::AndLeave);
    let mbar = MultiProgress::new();
    let factory = || mbar.add(ProgressBar::new_spinner().with_style(style.clone()));

    let pbar = factory().with_style(
        ProgressStyle::default_bar().template("[{elapsed_precise}] {msg}").on_finish(ProgressFinish::AndLeave),
    );
    pbar.set_message("Running...");

    // Parse core arguments and determine subcommand
    #[allow(clippy::type_complexity)]
    let (args, func): (&ArgMatches, Box<dyn FnOnce(&ArgMatches, CoreArgs) + Send>) = match app.subcommand() {
        Some(("rois", matches)) => (matches, Box::new(|matches, core| cli::roi::run(matches, core, factory))),
        Some(("loci", matches)) => (matches, Box::new(|matches, core| cli::loci::run(matches, core, factory))),
        _ => panic!("Subcommand is not specified."),
    };
    let core = cli::shared::args::CoreArgs::new(args, factory);

    // + 1 thread to render progress bar
    ThreadPoolBuilder::new().num_threads(core.threads + 1).build_global().expect(CREATE_THREAD_POOL_ERROR);
    rayon::scope(|s| {
        // Render progress bar in the additional thread
        s.spawn(|_| mbar.join().expect(RENDER_PROGRESS_ERROR));

        func(args, core);
        pbar.finish_with_message("Finished!");
    });
}
