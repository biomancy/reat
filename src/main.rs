use std::sync::{Arc, Mutex};

use clap::{crate_authors, crate_name, crate_version, App, AppSettings, ArgMatches};
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use rayon::ThreadPoolBuilder;

use rada::cli;
use rada::cli::shared::args::CoreArgs;

const CREATE_THREAD_POOL_ERROR: &str = "Failed to initialize thread pool";
const RENDER_PROGRESS_ERROR: &str = "Failed to render progress bar";

#[derive(Clone)]
struct PanicAwareProgressManager {
    mbar: Arc<MultiProgress>,
    lock: Arc<Mutex<Vec<ProgressBar>>>,
}

impl PanicAwareProgressManager {
    pub fn new() -> Self {
        let obj = Self { mbar: Arc::from(MultiProgress::new()), lock: Arc::new(Mutex::new(vec![])) };

        let hook = std::panic::take_hook();
        let clone = obj.clone();
        std::panic::set_hook(Box::new(move |info| {
            for x in clone.lock.lock().unwrap().iter() {
                x.abandon();
            }
            hook(info);
        }));
        obj
    }

    pub fn attach(&self, style: ProgressStyle) -> ProgressBar {
        let pbar = ProgressBar::new_spinner().with_style(style);
        self.mbar.add(pbar.clone());
        self.lock.lock().unwrap().push(pbar.clone());
        pbar
    }
}

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
                .args(cli::rois::args()),
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
    let masterbar = PanicAwareProgressManager::new();
    let factory = || masterbar.attach(style.clone());

    let pbar = masterbar.attach(
        ProgressStyle::default_spinner()
            .template("[{elapsed_precise}] {msg:.red.bold}")
            .on_finish(ProgressFinish::AndLeave),
    ); //.with_style(style);
    pbar.set_message("Running...");
    // Parse core arguments and determine subcommand
    #[allow(clippy::type_complexity)]
    let (args, func): (&ArgMatches, Box<dyn FnOnce(&ArgMatches, CoreArgs) + Send>) = match app.subcommand() {
        Some(("rois", matches)) => (matches, Box::new(|matches, core| cli::rois::run(matches, core, factory))),
        Some(("loci", matches)) => (matches, Box::new(|matches, core| cli::loci::run(matches, core, factory))),
        _ => panic!("Subcommand is not specified."),
    };
    let core = cli::shared::args::CoreArgs::new(args, factory);

    // + 1 thread to render progress bar
    ThreadPoolBuilder::new().num_threads(core.threads + 1).build_global().expect(CREATE_THREAD_POOL_ERROR);
    rayon::scope(|s| {
        // Render progress bar in the additional thread
        s.spawn(|_| {
            masterbar.mbar.join().expect(RENDER_PROGRESS_ERROR);
        });

        func(args, core);
        pbar.finish_with_message("Finished!");
    });
}
