use std::str::FromStr;

use crate::rada::stranding::deduct::StrandSpecificExperimentDesign;
use crate::rada::stranding::deduct::StrandSpecificExperimentDesign::*;

pub enum Stranding {
    Unstranded,
    Stranded(StrandSpecificExperimentDesign),
}

impl FromStr for Stranding {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "u" | "unstranded" => Ok(Stranding::Unstranded),
            "s" => Ok(Stranding::Stranded(Same)),
            "f" => Ok(Stranding::Stranded(Flip)),
            "s/f" => Ok(Stranding::Stranded(Same1Flip2)),
            "f/s" => Ok(Stranding::Stranded(Flip1Same2)),
            _ => Err(format!("Unknown strand: {}", s)),
        }
    }
}
