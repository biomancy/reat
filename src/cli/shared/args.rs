use clap::{Arg, ArgSettings};

use super::validate;

pub fn reqdefaults() -> Vec<ArgSettings> {
    vec![ArgSettings::Required, ArgSettings::TakesValue]
}

pub fn defaults() -> Vec<ArgSettings> {
    vec![ArgSettings::TakesValue]
}

pub mod core {
    use super::*;
    pub const INPUT: &str = "input";
    pub const REFERENCE: &str = "reference";
    pub const BINSIZE: &str = "binsize";
    pub const STRANDING: &str = "stranding";
    pub const THREADS: &str = "threads";
    pub const SAVETO: &str = "saveto";
    pub const NAME: &str = "name";

    pub const SECTION_NAME: &str = "Core";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(INPUT)
                .short('i')
                .long(INPUT)
                .settings(&reqdefaults())
                .multiple(true)
                .validator(validate::path)
                .long_about("Path to the input BAM file(s). May contain a space-separated list of files, in which case they are treated as technical replicates and pulled together."),
            Arg::new(REFERENCE)
                .short('r')
                .long(REFERENCE)
                .settings(&reqdefaults())
                .validator(validate::path)
                .long_about("Indexed fasta file with a reference genome assembly. Contig/chromosome names must match the names in the input BAM headers."),
            Arg::new(BINSIZE)
                .long(BINSIZE)
                .settings(&defaults())
                .validator(validate::numeric(1u32, 1_000_000u32))
                .default_value("64000")
                .long_about("Summarize the mismatches for each locus, providing each worker thread with genome bins(job share) of approximately the specified size (in base pairs)."),
            Arg::new(STRANDING)
                .short('s')
                .long(STRANDING)
                .settings(&reqdefaults())
                .validator(validate::stranding)
                .possible_values(&["u", "s", "f", "s/f", "f/s"])
                .long_about("Strand-specificity of the experiment, i.e. matching between the read strand and the transcript strand. Use \"u\" for unstranded experiments; other available options based on the RSeQC nomenclature(see infer_experiment.py docs): same:\"s\" (++,--), flip:\"f\" (+-,-+), same read1/flip read2:\"sf\" (1++,1--/2+-,2-+), flip read1/same read2:\"fs\" (1+-,1-+/2++,2--)."),
            Arg::new(NAME)
                .short('n')
                .long(NAME)
                .settings(&defaults())
                .default_value("NA")
                .long_about("Name of the run."),
            Arg::new(SAVETO)
                .short('o')
                .long(SAVETO)
                .settings(&defaults())
                .validator(validate::writable)
                .default_value("/dev/stdout")
                .long_about("Path to the output tsv file. By default, the results are printed to stdout."),
            Arg::new(THREADS)
                .short('t')
                .long(THREADS)
                .settings(&defaults())
                .validator(validate::numeric(1, usize::MAX))
                .default_value("1")
                .long_about("Maximum number of threads to spawn at once."),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod reads_filtering {
    use super::*;
    pub const MAPQ: &str = "mapq";
    pub const ALLOW_MAPQ_255: &str = "mapq-255";
    pub const INCLUDE_FLAGS: &str = "in-flags";
    pub const EXCLUDE_FLAGS: &str = "ex-flags";
    pub const PHREAD: &str = "phread";

    pub const SECTION_NAME: &str = "Reads filtering";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(MAPQ)
                .long(MAPQ)
                .settings(&defaults())
                .validator(validate::numeric(0u8, 254u8))
                .default_value("1")
                .long_about("Count only reads with mapq ≥ threshold. Note that reads with mapq = 255 are skipped by default(mapq is not available according to the SAM spec)."),
            Arg::new(ALLOW_MAPQ_255)
                .long(ALLOW_MAPQ_255)
                .settings(&defaults())
                .takes_value(false)
                .long_about("Count reads with mapq=255. Useful for aligners that do not fully conform to the SAM specification (e.g. STAR with default parameters)"),
            Arg::new(INCLUDE_FLAGS)
                .long(INCLUDE_FLAGS)
                .settings(&defaults())
                .validator(validate::numeric(0u16, 4095u16))
                .default_value("0")
                .long_about("Include only reads for which all the specified BAM flags are set. For example, a value of 3 will result in skipping reads that were not mapped in proper pairs. Use zero(0) to disable this filter."),
            Arg::new(EXCLUDE_FLAGS)
                .long(EXCLUDE_FLAGS)
                .settings(&defaults())
                .validator(validate::numeric(0u16, 4095u16))
                .default_value("2820")
                .long_about("Exclude reads for which any of the specified BAM flags are set. For example, a value of 2820 will result in skipping unmapped reads, supplementary and secondary alignments, reads that fail platform/vendor quality checks. Use zero(0) to disable this filter."),
            Arg::new(PHREAD)
                .long(PHREAD)
                .settings(&defaults())
                .validator(validate::numeric(0u8, 255u8))
                .default_value("20")
                .long_about("Count only bases with phread ≥ threshold. For reference, Phread is defined as -10 log₁₀[error probability], so Phread = 20 means 1 error in 100 base calls."),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod autoref {
    use super::*;
    pub const MIN_COVERAGE: &str = "ref-min-cov";
    pub const MIN_FREQ: &str = "ref-min-freq";
    pub const HYPEREDITING: &str = "hyperedit";

    pub const SECTION_NAME: &str = "Autoref";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(MIN_COVERAGE)
                .long(MIN_COVERAGE)
                .settings(&defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("20")
                .long_about("Automatically correct reference sequence for loci with coverage ≥ the threshold. In short, there is no reason to use the assembly nucleotide \"T \" if we have sequenced 100% \"A \". This heuristic is especially useful in regions of low complexity(or simple repeats), where such SNPs can affect the editing estimation."),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .settings(&defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.95")
                .long_about("Automatically correct reference sequence for loci with the most common nucleotide frequency ≥ cutoff"),
            Arg::new(HYPEREDITING)
                .long(HYPEREDITING)
                .settings(&defaults())
                .takes_value(false)
                .long_about("Turn on the \"hyperediting\" mode, i.e. do not correct(replace) A with G and T with C. This will ensure that potentially hyper-editable sites are not accidentally lost.")
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod stranding {
    use super::*;
    pub const MIN_MISMATCHES: &str = "str-min-mismatches";
    pub const MIN_FREQ: &str = "str-min-freq";
    pub const ANNOTATION: &str = "annotation";

    pub const SECTION_NAME: &str = "Stranding";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(ANNOTATION)
                .long(ANNOTATION)
                .settings(&defaults())
                .validator(validate::path)
                .long_about("Genome annotation in the GFF3 format. Genomic features (exons and genes) are used to inference loci/ROI strand based on the most likely direction of transcription (see the GitHub documentation for details). It is recommended to provide genome annotation for unstranded libraries, otherwise stranding will be highly inaccurate."),
            Arg::new(MIN_MISMATCHES)
                .long(MIN_MISMATCHES)
                .settings(&defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("50")
                .long_about("Automatically predict strand based on the observed A->I editing for locus/ROI with A->G mismatches >= threshold. It is a fallback strand prediction heuristic, used only for the unstranded libraries. Not relevant for organisms without active ADAR-like enzymes."),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .settings(&defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.05")
                .long_about("Automatically predict strand based on the observed A->I editing for locus/ROI with A->G freq >= threshold (freq = ∑ A->G / (∑ A->G + ∑ A->A))"),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod output_filtering {
    use super::*;

    pub const MIN_MISMATCHES: &str = "out-min-mismatches";
    pub const MIN_FREQ: &str = "out-min-freq";

    pub const SECTION_NAME: &str = "Output filtering";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(MIN_MISMATCHES)
                .long(MIN_MISMATCHES)
                .settings(&defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("5")
                .long_about("Output only ROI having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". For \"N\" reference all nucleotides are considered as mismatches. This is a deliberate choice to allow a subsequent user to work through / filter such records."),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .settings(&defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.01")
                .long_about("Output only ROI having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)"),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub fn all<'a>() -> Vec<Arg<'a>> {
    core::args()
        .into_iter()
        .chain(reads_filtering::args().into_iter())
        .chain(autoref::args().into_iter())
        .chain(stranding::args().into_iter())
        .collect()
}

use std::path::PathBuf;

use clap::ArgMatches;

use indicatif::ProgressBar;
use rust_htslib::bam::Record;

use crate::cli::shared::stranding::Stranding;
use crate::core::filtering::reads::{ReadsFilterByFlags, ReadsFilterByQuality, SequentialReadsFilter};
use crate::core::refnuc::RefNucPredByHeurisitc;

use super::parse;
use std::fs::File;
use std::io::BufWriter;

pub struct CoreArgs {
    pub name: String,
    pub threads: usize,
    pub bamfiles: Vec<PathBuf>,
    pub reference: PathBuf,
    pub refnucpred: RefNucPredByHeurisitc,
    pub readfilter: SequentialReadsFilter<Record, ReadsFilterByQuality, ReadsFilterByFlags>,
    pub stranding: Stranding,
    pub saveto: BufWriter<File>,
}

impl CoreArgs {
    pub fn new(args: &ArgMatches, factory: impl Fn() -> ProgressBar) -> Self {
        Self {
            name: parse::name(factory(), args),
            threads: parse::threads(factory(), args),
            bamfiles: parse::bamfiles(factory(), args),
            reference: parse::reference(factory(), args),
            refnucpred: parse::refnucpred(factory(), args),
            readfilter: parse::readfilter(factory(), args),
            stranding: parse::stranding(factory(), args),
            saveto: parse::saveto(factory(), args),
        }
    }
}
