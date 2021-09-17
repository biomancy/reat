


use clap::ArgMatches;
use clap::{Arg};

use indicatif::{ProgressBar};


use crate::cli::shared;
use crate::cli::shared::args::{defaults};

use crate::cli::shared::validate;

use crate::core::filtering::summary::SummaryFilterByMismatches;

use crate::core::stranding::predict::SequentialStrandPredictor;
use crate::core::workload::ROIWorkload;

use super::parse;

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
                .long_about("Output only Loci having total number of mismatches ≥ threshold. Mismatches are counted jointly, i.e. for the \"A\" reference we have \"C\" + \"G\" + \"T\". For \"N\" reference all nucleotides are considered as mismatches. This is a deliberate choice to allow a subsequent user to work through / filter such records."),
            Arg::new(MIN_FREQ)
                .long(MIN_FREQ)
                .settings(&defaults())
                .validator(validate::numeric(0f32, 1f32))
                .default_value("0.01")
                .long_about("Output only Loci having total mismatches frequency ≥ threshold (freq = ∑ mismatches / coverage)"),
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
    pub outfilter: SummaryFilterByMismatches,
    pub strandpred: SequentialStrandPredictor,
}

impl LociArgs {
    pub fn new(core: &shared::args::CoreArgs, args: &ArgMatches, factory: &impl Fn() -> ProgressBar) -> Self {
        let outfilter =
            shared::parse::outfilter(factory(), output_filtering::MIN_MISMATCHES, output_filtering::MIN_FREQ, args);

        let mut strandpred: Option<SequentialStrandPredictor> = Default::default();
        let mut workload: Option<Vec<ROIWorkload>> = Default::default();
        let mut maxsize: Option<u32> = Default::default();

        let (pbarw, pbars) = (factory(), factory());
        rayon::scope(|s| {
            s.spawn(|_| {
                let (w, m) = parse::work(pbarw, &core.bamfiles, args);
                workload = Some(w);
                maxsize = Some(m)
            });
            s.spawn(|_| strandpred = Some(shared::parse::strandpred(pbars, args)));
        });

        Self { workload: workload.unwrap(), maxwsize: maxsize.unwrap(), outfilter, strandpred: strandpred.unwrap() }
    }
}
