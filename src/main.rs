use clap::{crate_authors, crate_name, crate_version, App, AppSettings};

mod cli;
mod rada;

fn main() {
    let matches = App::new(crate_name!())
        .author(crate_authors!("\n"))
        .version(crate_version!())
        .max_term_width(120)
        .setting(AppSettings::DeriveDisplayOrder)
        .args(cli::args::all())
        .get_matches();

    cli::App::new(&matches).run();
}
