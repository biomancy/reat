use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use clap::ArgMatches;
use clap::{Arg, ArgSettings};
use derive_getters::Dissolve;
use indicatif::{MultiProgress, ProgressBar};
use rust_htslib::bam::Record;

use crate::cli::shared;
use crate::cli::shared::args::{defaults, reqdefaults};
use crate::cli::shared::stranding::Stranding;
use crate::cli::shared::validate;
use crate::core::filtering::reads::{ReadsFilterByFlags, ReadsFilterByQuality, SequentialReadsFilter};
use crate::core::filtering::summary::SummaryFilterByMismatches;
use crate::core::refnuc::RefNucPredByHeurisitc;
use crate::core::stranding::predict::SequentialStrandPredictor;
use crate::core::workload::ROIWorkload;

use super::parse;

pub mod stats {
    use super::*;

    pub const EDITING_INDEX: &str = "ei";

    pub const SECTION_NAME: &str = "Stats";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(EDITING_INDEX)
                .long(EDITING_INDEX)
                .settings(&defaults())
                .validator(validate::writable)
                .long_about("File for saving Editing Indexes (EI). If the file already exists, the EI results for the current experiments will be appended to it.")
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod special {
    use super::*;

    pub const ROI: &str = "roi";

    pub const SECTION_NAME: &str = "Special information";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(ROI)
                .long(ROI)
                .settings(&reqdefaults())
                .validator(validate::path)
                .long_about("Path to a BED-like file with 4 columns (chr, start, end, name) in which the target regions of interest (ROI) are declared. If specified, total mismatches will be reported for each region of interest rather than loci."),
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
    shared::args::all()
        .into_iter()
        .chain(stats::args())
        .chain(special::args())
        .chain(output_filtering::args())
        .collect()
}

pub struct ROIArgs {
    pub workload: Vec<ROIWorkload>,
    pub maxwsize: u32,
    pub outfilter: SummaryFilterByMismatches,
    pub strandpred: SequentialStrandPredictor,
    pub ei: Option<BufWriter<File>>,
}

impl ROIArgs {
    pub fn new(_: &shared::args::CoreArgs, args: &ArgMatches, factory: &impl Fn() -> ProgressBar) -> Self {
        let outfilter =
            shared::parse::outfilter(factory(), output_filtering::MIN_MISMATCHES, output_filtering::MIN_FREQ, args);
        let ei = parse::editing_index(factory(), args);

        let mut strandpred: Option<SequentialStrandPredictor> = Default::default();
        let mut workload: Option<Vec<ROIWorkload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();

        let (pbarw, pbars) = (factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::work(pbarw, args);
                workload = Some(w);
                maxsize = Some(m)
            });
            s.spawn(|_| strandpred = Some(shared::parse::strandpred(pbars, args)));
        });

        Self { workload: workload.unwrap(), maxwsize: maxsize.unwrap(), outfilter, strandpred: strandpred.unwrap(), ei }
    }
}
