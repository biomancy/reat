use clap::App;

use rada::cli;

mod paths;

fn tmp() {
    // let saveto = tempfile::NamedTempFile::new().unwrap();
    // let ei = tempfile::NamedTempFile::new().unwrap();
    //
    // let args = format!(
    //     "rada --input {input} --roi {roi} --reference {fasta} --annotation {gtf} \
    //     --saveto {saveto} --stranding u --threads 1 --ei {ei}",
    //     input = "",
    //     roi = "",
    //     fasta = "",
    //     gtf = paths::GRCh38::gtf().display(),
    //     saveto = saveto.path().display(),
    //     ei = ei.path().display()
    // );
    //
    // let matches = App::new("").args(cli::args::all()).get_matches_from(args);
    // cli::App::new(&matches).run();
    //
    // saveto.close().unwrap();
}
