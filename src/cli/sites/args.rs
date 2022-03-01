use std::collections::{HashMap, HashSet};

use clap::Arg;
use clap::ArgMatches;
use indicatif::ProgressBar;

use super::parse;
use crate::cli::shared;
use crate::cli::shared::args::defaults;
use crate::cli::shared::validate;
use crate::core::mismatches::prefilters;
use crate::core::workload::ROIWorkload;

pub mod output_filtering {
    use super::*;

    pub const MIN_MISMATCHES: &str = "out-min-mismatches";
    pub const MIN_FREQ: &str = "out-min-freq";
    pub const MIN_COVERAGE: &str = "out-min-cov";
    pub const FORCE_LIST: &str = "force";

    pub const SECTION_NAME: &str = "Output hooks";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args = vec![
            Arg::new(MIN_COVERAGE)
                .long(MIN_COVERAGE)
                .setting(defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("10")
                .long_help("Output only site covered by at least X unique filters(after filters/bases hooks)"),
            Arg::new(MIN_MISMATCHES)
                .long(MIN_MISMATCHES)
                .setting(defaults())
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("3")
                .long_help(
                    "Output only site having total number of mismatches ≥ threshold. \
                    Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". \
                    For \"N\" reference all nucleotides are considered as mismatches. \
                    This is a deliberate choice to allow a subsequent user to work through / filter such records.",
                ),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .setting(defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.01")
                .long_help(
                    "Output only site having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)",
                ),
            Arg::new(FORCE_LIST)
                .long(FORCE_LIST)
                .setting(defaults())
                .validator(validate::path)
                .long_help(
                    "Force the output of sites overlapping regions in a given BED file, even if they do not pass other filters"
                ),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub fn all<'a>() -> Vec<Arg<'a>> {
    shared::args::all().into_iter().chain(output_filtering::args()).collect()
}

pub struct LociArgs {
    pub workload: Vec<ROIWorkload>,
    pub maxwsize: u32,
    pub filter: prefilters::ByMismatches,
    // pub stranding: IntervalStrandingEngine<BorrowedIntervalMismatches<'a>>,
    // pub forcelist: Option<HashMap<String, HashSet<u64>>>,
}

impl LociArgs {
    pub fn new(core: &shared::args::CoreArgs, args: &ArgMatches, factory: &impl Fn() -> ProgressBar) -> Self {
        let filter = parse::outfilter(
            factory(),
            output_filtering::MIN_MISMATCHES,
            output_filtering::MIN_FREQ,
            output_filtering::MIN_COVERAGE,
            args,
        );

        // let mut stranding: Option<IntervalStrandingEngine<BorrowedIntervalMismatches<'_>>> = Default::default();
        let mut workload: Option<Vec<ROIWorkload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();
        // let mut forcelist: Option<HashMap<String, HashSet<u64>>> = Default::default();

        let (pbarw, pbars, pbarf) = (factory(), factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::work(pbarw, &core.bamfiles, args);
                workload = Some(w);
                maxsize = Some(m)
            });
            // s.spawn(|_| stranding = Some(shared::parse::stranding(pbars, args)));
            // s.spawn(|_| forcelist = parse::forcein(pbarf, args))
        });

        Self {
            workload: workload.unwrap(),
            maxwsize: maxsize.unwrap(),
            filter,
            // stranding: stranding.unwrap(),
            // forcelist,
        }
    }
}
