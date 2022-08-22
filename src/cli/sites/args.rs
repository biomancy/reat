use clap::Arg;
use clap::ArgMatches;
use indicatif::ProgressBar;

use crate::cli::shared;
use crate::cli::shared::validate;
use crate::core::mismatches::prefilters;
use crate::core::mismatches::prefilters::retain::RetainSitesFromIntervals;
use crate::core::mismatches::site::SiteMismatchesVec;
use crate::core::stranding::predict::REATStrandingEngine;
use crate::core::workload::SiteWorkload;

use super::parse;

pub mod output_filtering {
    use super::*;

    pub const MIN_MISMATCHES: &str = "out-min-mismatches";
    pub const MIN_FREQ: &str = "out-min-freq";
    pub const MIN_COVERAGE: &str = "out-min-cov";
    pub const FORCE_LIST: &str = "force";
    pub const REGIONS: &str = "region";

    pub const SECTION_NAME: &str = "Output hooks";

    pub fn args<'a>() -> Vec<Arg<'a>> {
        let args =
            vec![
            Arg::new(MIN_COVERAGE)
                .long(MIN_COVERAGE)
                .takes_value(true)
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("10")
                .long_help("Output only site covered by at least X unique reads"),
            Arg::new(MIN_MISMATCHES)
                .long(MIN_MISMATCHES)
                .takes_value(true)
                .validator(validate::numeric(0u32, u32::MAX))
                .default_value("3")
                .long_help(
                    "Output only sites with total number of mismatches ≥ threshold. \
                    Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". \
                    For \"N\" reference all nucleotides are considered as mismatches. \
                    This is a deliberate choice to allow a subsequent user to work through / filter such records.",
                ),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .takes_value(true)
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.01")
                .long_help(
                    "Output only sites with total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)",
                ),
            Arg::new(FORCE_LIST).long(FORCE_LIST).takes_value(true).validator(validate::path).long_help(
                "Force the output of sites located in a given BED file (even if they do not pass other filters).",
            ),
            Arg::new(REGIONS).long(REGIONS).takes_value(true).validator(validate::path).long_help(
                "Process only sites overlapping the given BED file.",
            ),
        ];
        args.into_iter().map(|x| x.help_heading(Some(SECTION_NAME))).collect()
    }
}

pub fn all<'a>() -> Vec<Arg<'a>> {
    shared::args::all().into_iter().chain(output_filtering::args()).collect()
}

pub struct SiteArgs {
    pub workload: Vec<SiteWorkload>,
    pub maxwsize: usize,
    pub prefilter: prefilters::ByMismatches,
    pub stranding: REATStrandingEngine<SiteMismatchesVec>,
    pub retain: Option<RetainSitesFromIntervals>,
}

impl SiteArgs {
    pub fn new(core: &mut shared::args::CoreArgs, args: &ArgMatches, factory: &impl Fn() -> ProgressBar) -> Self {
        let filter = shared::parse::outfilter(
            factory(),
            output_filtering::MIN_MISMATCHES,
            output_filtering::MIN_FREQ,
            output_filtering::MIN_COVERAGE,
            args,
        );

        let mut stranding = REATStrandingEngine::new();
        let mut workload: Option<Vec<SiteWorkload>> = Default::default();
        let mut maxsize: Option<usize> = Default::default();
        let mut retain: Option<RetainSitesFromIntervals> = Default::default();

        let (pbarw, pbars, pbarf) = (factory(), factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::work(pbarw, &core.bamfiles, core.excluded.take(), args);
                workload = Some(w);
                maxsize = Some(m)
            });
            s.spawn(|_| {
                stranding = shared::parse::strandpred(pbars, args);
            });
            s.spawn(|_| retain = parse::retain(pbarf, args));
        });

        Self { workload: workload.unwrap(), maxwsize: maxsize.unwrap(), prefilter: filter, stranding, retain }
    }
}
