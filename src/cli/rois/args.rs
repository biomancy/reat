use std::fs::File;
use std::io::BufWriter;

use clap::Arg;
use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::cli::shared::args::{defaults, reqdefaults};
use crate::cli::shared::validate;
use crate::core::hooks::filters::ByMismatches;
use crate::core::stranding::predict::StrandingEngine;
use crate::core::workload::ROIWorkload;

use super::parse;
use crate::cli::shared::stranding::Stranding;
use crate::core::hooks::stats::ROIEditingIndex;
use crate::core::mismatches::IntermediateMismatches;
use crate::core::stranding::predict::algo::{StrandByAtoIEditing, StrandByGenomicAnnotation};
use crate::core::stranding::predict::engines::ROIStrandingEngine;

pub mod stats {
    use super::*;

    pub const EDITING_INDEX: &str = "ei";

    pub const SECTION_NAME: &str = "Stats";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![Arg::new(EDITING_INDEX)
            .long(EDITING_INDEX)
            .setting(defaults())
            .validator(validate::writable)
            .long_help(
                "File for saving Editing Indexes (EI). \
                If the file already exists, EI for the current experiments will be appended to it",
            )];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod special {
    use super::*;

    pub const ROI: &str = "rois";

    pub const SECTION_NAME: &str = "Special information";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![Arg::new(ROI).long(ROI).setting(reqdefaults()).validator(validate::path).long_help(
            "Path to a BED file with regions of interest(ROIS) \
            with at least 4 first BED columns(chr, start, end, name)",
        )];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub mod output_filtering {
    use super::*;

    pub const MIN_MISMATCHES: &str = "out-min-mismatches";
    pub const MIN_FREQ: &str = "out-min-freq";
    pub const MIN_COVERAGE: &str = "out-min-cov";

    pub const SECTION_NAME: &str = "Output hooks";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(MIN_COVERAGE)
                .long(MIN_COVERAGE)
                .setting(defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("10")
                .long_help("Output only ROIs covered by at least X unique filters(after filters/bases hooks)"),
            Arg::new(MIN_MISMATCHES)
                .long(MIN_MISMATCHES)
                .setting(defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("5")
                .long_help(
                    "Output only ROI having total number of mismatches ≥ threshold. \
                    Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". \
                    For \"N\" reference all nucleotides are considered as mismatches. \
                    This is a deliberate choice to allow a subsequent user to work through / filter such records",
                ),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .setting(defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.01")
                .long_help(
                    "Output only ROI having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)",
                ),
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
    pub filter: ByMismatches,
    pub ei: Option<BufWriter<File>>,
}

impl ROIArgs {
    pub fn new(core: &shared::args::CoreArgs, args: &ArgMatches, factory: &impl Fn() -> ProgressBar) -> Self {
        let filter = shared::parse::outfilter(
            factory(),
            output_filtering::MIN_MISMATCHES,
            output_filtering::MIN_FREQ,
            output_filtering::MIN_COVERAGE,
            args,
        );
        let ei = parse::editing_index(factory(), args);

        // let mut stranding: Option<Stranding> = Default::default();
        let mut workload: Option<Vec<ROIWorkload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();

        let (pbarw, pbars) = (factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::work(pbarw, args, core.excluded.clone());
                workload = Some(w);
                maxsize = Some(m)
            });
            // s.spawn(|_| strandpred = Some(shared::parse::stranding(pbars, args)));
        });

        Self { workload: workload.unwrap(), maxwsize: maxsize.unwrap(), filter, ei }
    }
}
