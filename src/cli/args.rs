use clap::{Arg, ArgSettings};

use super::validate;

// Core
pub const INPUT: &'static str = "input";
pub const REFERENCE: &'static str = "reference";
pub const ROI: &'static str = "roi";
pub const BINSIZE: &'static str = "binsize";
pub const STRANDING: &'static str = "stranding";
pub const THREADS: &'static str = "threads";
pub const SAVETO: &'static str = "saveto";

// Filtering
pub const MAPQ: &'static str = "mapq";
pub const PHREAD: &'static str = "phread";
pub const OUT_MIN_MISMATCHES: &'static str = "out-min-mismatches";
pub const OUT_MIN_FREQ: &'static str = "out-min-freq";

// Autoref
pub const AUTOREF_MIN_COVERAGE: &'static str = "ref-min-cov";
pub const AUTOREF_MIN_FREQ: &'static str = "ref-min-freq";

// Stranding
pub const STRANDING_MIN_MISMATCHES: &'static str = "str-min-mismatches";
pub const STRANDING_MIN_FREQ: &'static str = "str-min-freq";
pub const ANNOTATION: &'static str = "annotation";

fn reqdefaults() -> Vec<ArgSettings> {
    vec![ArgSettings::Required, ArgSettings::TakesValue, ArgSettings::AllowHyphenValues]
}

fn defaults() -> Vec<ArgSettings> {
    vec![ArgSettings::TakesValue, ArgSettings::AllowHyphenValues]
}

pub fn core<'a>() -> Vec<Arg<'a>> {
    let args = vec![
        Arg::new(INPUT)
            .short('i')
            .long(INPUT)
            .settings(&reqdefaults())
            .multiple(true)
            .validator(validate::path)
            .long_about("Path to the input BAM file(s). May contain a comma-separated list of files, in which case they are treated as technical replicates and pulled together."),
        Arg::new(REFERENCE)
            .short('r')
            .long(REFERENCE)
            .settings(&reqdefaults())
            .validator(validate::path)
            .long_about("Indexed fasta file with a reference genome assembly. Contig/chromosome names must match the names in the input BAM headers."),
        Arg::new(ROI)
            .long(ROI)
            .settings(&reqdefaults())
            .validator(validate::path)
            .long_about("Path to a BED-like file with 4 columns (chr, start, end, name) in which the target regions of interest (ROI) are declared. If specified, cumulative mismatches will be reported for each region of interest rather than loci.")
            .conflicts_with(BINSIZE),
        Arg::new(BINSIZE)
            .long(BINSIZE)
            .settings(&reqdefaults())
            .validator(validate::numeric(1u32, 1_000_000u32))
            .long_about("Summarize the mismatches for each locus, providing each worker thread with genome bins(job share) of approximately the specified size (in base pairs)."),
        Arg::new(STRANDING)
            .short('s')
            .long(STRANDING)
            .settings(&reqdefaults())
            .validator(validate::stranding)
            .long_about("Strand-specificity of the experiment, i.e. matching the read strand and the transcript strand. Use \"u\" for unstranded experiments; other available options based on the RSeQC nomenclature(see infer_experiment.py docs): \"s\" (++,--), \"f\" (+-,-+), \"sf\" (1++,1--,2+-,2-+), \"fs\" (1+-,1-+,2++,2--)."),
        Arg::new(SAVETO)
            .short('o')
            .long(SAVETO)
            .settings(&reqdefaults())
            .validator(validate::writable)
            .long_about("Path to the output tsv file. May end with \".gz\", in which case the stream is automatically gzipped."),
        Arg::new(THREADS)
            .short('t')
            .long(THREADS)
            .settings(&defaults())
            .default_value("1")
            .long_about("Maximum number of threads to spawn at once."),
    ];
    args.into_iter().map(|x| x.help_heading(Some("Core"))).collect()
}

pub fn filtering<'a>() -> Vec<Arg<'a>> {
    let args = vec![
        Arg::new(MAPQ)
            .long(MAPQ)
            .settings(&defaults())
            .validator(validate::numeric(0u8, 254u8))
            .default_value("1")
            .long_about("Count only reads with mapq ≥ threshold. Note that reads with mapq = 255 are always skipped(mapq is not available according to the SAM spec)."),
        Arg::new(PHREAD)
            .long(PHREAD)
            .settings(&defaults())
            .validator(validate::numeric(0u8, 255u8))
            .default_value("20")
            .long_about("Count only bases with phread ≥ threshold. For reference, Phread is defined as -10 log₁₀[error probability], so Phread = 20 means 1 error in 100 base calls."),
        Arg::new(OUT_MIN_MISMATCHES)
            .long(OUT_MIN_MISMATCHES)
            .settings(&defaults())
            .validator(validate::numeric(0u32, u32::MAX))
            .default_value("20")
            .long_about("Output only loci/ROI having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". For \"N\" reference all nucleotides are considered as mismatches. This is a deliberate choice to allow a subsequent user to work through / filter such records."),
        Arg::new(OUT_MIN_FREQ)
            .long(OUT_MIN_FREQ)
            .settings(&defaults())
            .validator(validate::numeric(0f32, 1f32))
            .default_value("0.05")
            .long_about("Output only loci/ROI having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)"),
    ];
    args.into_iter().map(|x| x.help_heading(Some("Filtering"))).collect()
}

pub fn autoref<'a>() -> Vec<Arg<'a>> {
    let args = vec![
        Arg::new(AUTOREF_MIN_COVERAGE)
            .long(AUTOREF_MIN_COVERAGE)
            .settings(&defaults())
            .validator(validate::numeric(0u32, u32::MAX))
            .default_value("20")
            .long_about("Automatically correct reference sequence for loci with coverage ≥ the threshold. In short, there is no reason to use the assembly nucleotide \"T \" if we have sequenced 100% \"A \". This heuristic is especially useful in regions of low complexity(or simple repeats), where such SNPs can affect the editing estimation."),
        Arg::new(AUTOREF_MIN_FREQ)
            .long(AUTOREF_MIN_FREQ)
            .settings(&defaults())
            .validator(validate::numeric(0f32, 1f32))
            .default_value("0.95")
            .about("Automatically correct reference sequence for loci with the most common nucleotide frequency ≥ cutoff"),
    ];
    args.into_iter().map(|x| x.help_heading(Some("Autoref"))).collect()
}

pub fn stranding<'a>() -> Vec<Arg<'a>> {
    let args = vec![
        Arg::new(ANNOTATION)
            .long(ANNOTATION)
            .settings(&defaults())
            .validator(validate::path)
            .long_about("Genome annotation in the GFF3 format. Genomic features (exons and genes) are used to inference loci/ROI strand based on the most likely direction of transcription (see the GitHub documentation for details). It is highly recommended to provide genome annotation for unstranded libraries, otherwise stranding will be highly inaccurate."),
        Arg::new(STRANDING_MIN_MISMATCHES)
            .long(STRANDING_MIN_MISMATCHES)
            .settings(&defaults())
            .validator(validate::numeric(0u32, u32::MAX))
            .default_value("50")
            .long_about("Automatically predict strand based on the observed A->I editing for locus/ROI with A->G mismatches >= threshold. It is a fallback strand prediction heuristic, used only for the unstranded libraries. Not relevant for organisms without active ADAR-like enzymes."),
        Arg::new(STRANDING_MIN_FREQ)
            .long(STRANDING_MIN_FREQ)
            .settings(&defaults())
            .validator(validate::numeric(0f32, 1f32))
            .default_value("0.05")
            .about("Automatically predict strand based on the observed A->I editing for locus/ROI with A->G frequency >= threshold (freq = ∑ A->G / (∑ A->G + ∑ A->A))"),
    ];
    args.into_iter().map(|x| x.help_heading(Some("Stranding"))).collect()
}

pub fn all<'a>() -> Vec<Arg<'a>> {
    core()
        .into_iter()
        .chain(filtering().into_iter())
        .chain(autoref().into_iter())
        .chain(stranding().into_iter())
        .collect()
}
