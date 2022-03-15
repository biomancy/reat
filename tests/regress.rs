// use clap::App;

// use rada::cli;

use std::fs::File;
use std::path::Path;

use clap::App;
use file_diff::diff_files;
use indicatif::{MultiProgress, ProgressBar};
use rayon::ThreadPoolBuilder;
use tempfile::NamedTempFile;

use reat::cli;

mod paths;

const TMP_CREATE_ERROR: &'static str = "Failed to create temporary file";
const TMP_DELETE_ERROR: &'static str = "Failed to delete temporary file";
const THREAD_POOL_ERROR: &'static str = "Failed to initialize thread pool";

#[allow(non_camel_case_types)]
enum SubCommand {
    sites,
    rois,
}

fn run(args: &[&str], launch: SubCommand) {
    let masterbar = MultiProgress::new();
    let factory = || masterbar.add(ProgressBar::hidden());

    let app = match launch {
        SubCommand::sites => cli::sites::args(),
        SubCommand::rois => cli::rois::args(),
    };

    let app = App::new("test").args(app);
    let args = app.get_matches_from(args);

    let core = cli::shared::args::CoreArgs::new(&args, factory);
    let pool = ThreadPoolBuilder::new().num_threads(core.threads).build().expect(THREAD_POOL_ERROR);
    pool.scope(|_| match launch {
        SubCommand::sites => cli::sites::run(&args, core, factory),
        SubCommand::rois => cli::rois::run(&args, core, factory),
    });
    masterbar.join_and_clear().expect("Failed to join pbars. Leak?");
}

fn same(first: &Path, second: &Path) -> bool {
    let mut first = match File::open(first) {
        Ok(f) => f,
        Err(e) => panic!("{}", e),
    };
    let mut second = match File::open(second) {
        Ok(f) => f,
        Err(e) => panic!("{}", e),
    };
    diff_files(&mut first, &mut second)
}

mod site {
    use super::*;

    #[test]
    fn trimming() {
        // ../target/release/reat site --input resources/bam/SRX6966474.bam \
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "f/s" \
        // -n Test -5 10 --trim3 2 -o resources/expected/site/trimmed.csv
        let expected = paths::expected::LOCI.join("trimmed.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "--input", &paths::bam::EXAMPLE, "--reference", &paths::GRCh38::FASTA, "-s", "f/s",
            "-n", "Test", "-5", "10", "--trim3", "2", "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::sites);

        assert!(same(tmp.path(), &expected));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn deducted_strand() {
        // ../target/release/reat site --input resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        // -s "f/s" -n Test --out-min-cov 20 -o resources/expected/site/deducted.csv
        let expected = paths::expected::LOCI.join("deducted.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "--input", &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA, "-s", "f/s",
            "-n", "Test", "--out-min-cov", "20", "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::sites);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn predicted_strand() {
        // ../target/release/reat site --input resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        // -s u --hyperedit --ref-min-cov 30 --annotation resources/GRCh38/Homo_sapiens.GRCh38.104.gff3.gz
        // --str-min-freq 0.01 --str-min-mismatches 5 -o resources/expected/site/predicted.csv
        let expected = paths::expected::LOCI.join("predicted.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, "--reference", &paths::GRCh38::FASTA, "--stranding", "u",
            "--hyperedit", "--ref-min-cov", "30", "--annotation", &paths::GRCh38::GFF3, "--str-min-freq", "0.01",
            "--str-min-mismatches", "5", "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::sites);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn multiple_files() {
        // ../target/release/reat site --input resources/bam/SRX6966474.bam resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s f -t 12 --in-flags 67
        // --out-min-freq 0.1 -o resources/expected/site/doubled.csv
        let expected = paths::expected::LOCI.join("doubled.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, &paths::bam::EXAMPLE, "--reference", &paths::GRCh38::FASTA,
            "--stranding", "f", "-t", "2", "--in-flags", "67", "--out-min-freq", "0.1",
            "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::sites);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }
}

mod rois {
    use std::fs;

    use super::*;

    #[test]
    fn trimming() {
        // ../target/release/reat roi -i resources/bam/SRX6966474.bam \
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "f/s" --out-min-cov 30 -n Test \
        // --in-flags=3 --trim5=2 -3 10 --rois resources/GRCh38/Alu.bed.gz --hyperedit \
        // -o resources/expected/rois/trimmed.csv
        let expected = paths::expected::ROIS.join("trimmed.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA, "-s", "f/s", "--out-min-cov", "30",
            "-n", "Test", "--in-flags", "3", "--trim5", "2", "-3", "10", "--rois", &paths::GRCh38::ALU,
            "--hyperedit", "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::rois);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn deducted_strand() {
        // ../target/release/reat roi -i resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "f/s" --out-min-cov 20 -n Test
        // --in-flags=3 --rois resources/GRCh38/Alu.bed.gz --hyperedit -o resources/expected/rois/deducted.csv
        let expected = paths::expected::ROIS.join("deducted.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA, "-s", "f/s", "--out-min-cov", "20",
            "-n", "Test", "--in-flags", "3", "--rois", &paths::GRCh38::ALU, "--hyperedit",
            "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::rois);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn predicted_strand() {
        // ../target/release/reat roi -i resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "u"
        // --rois resources/GRCh38/Alu.bed.gz --ref-min-cov=30
        // --annotation resources/GRCh38/Homo_sapiens.GRCh38.104.gff3.gz --str-min-freq 0.01
        // --str-min-mismatches 5 -o resources/expected/rois/predicted.csv
        let expected = paths::expected::ROIS.join("predicted.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA, "--stranding", "u",
            "--rois", &paths::GRCh38::ALU, "--ref-min-cov", "30",
            "--annotation", &paths::GRCh38::GFF3, "--str-min-freq", "0.01",
            "--str-min-mismatches", "5", "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::rois);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn multiple_files() {
        // ../target/release/reat roi -i resources/bam/SRX6966474.bam resources/bam/SRX6966474.bam
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "s/f" --ref-min-freq 0.8
        // --rois resources/GRCh38/Alu.bed.gz -o resources/expected/rois/doubled.csv -t 2
        let expected = paths::expected::ROIS.join("doubled.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        #[rustfmt::skip]
            let args = [
            "test", "-i", &paths::bam::EXAMPLE, &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA,
            "--stranding", "s/f", "--rois", &paths::GRCh38::ALU, "--ref-min-freq", "0.8", "-t", "2",
            "-o", tmp.path().to_str().unwrap(),
        ];
        run(&args, SubCommand::rois);

        assert!(same(tmp.path(), &expected.as_path()));
        tmp.close().expect(TMP_DELETE_ERROR);
    }

    #[test]
    fn ei() {
        // for name in "Test 1" ".";
        // do
        //  ../target/release/reat roi -i resources/bam/SRX6966474.bam \
        // -r resources/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -s "f/s" -n "$name" --in-flags=3 --rois \
        // resources/GRCh38/Alu.bed.gz -o /dev/null --ei resources/expected/rois/ei.csv; \
        // done

        let expected = paths::expected::ROIS.join("ei.csv");
        assert!(expected.is_file());

        let tmp = NamedTempFile::new().expect(TMP_CREATE_ERROR);
        let ei = tmp.path().to_owned();

        // Header will be printed only if an EI file doesn't exist
        tmp.close().expect(TMP_DELETE_ERROR);

        for name in ["Test 1", "."] {
            #[rustfmt::skip]
                let args = [
                "test", "-i", &paths::bam::EXAMPLE, "-r", &paths::GRCh38::FASTA, "-s", "f/s",
                "-n", name, "--in-flags", "3", "--rois", &paths::GRCh38::ALU,
                "-o", "/dev/null", "--ei", ei.to_str().unwrap()
            ];
            run(&args, SubCommand::rois);
        }

        assert!(same(ei.as_path(), &expected.as_path()), "{} vs {}", ei.display(), expected.display());
        fs::remove_file(ei).expect(TMP_DELETE_ERROR);
    }
}
