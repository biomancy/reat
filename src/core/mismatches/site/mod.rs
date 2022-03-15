use std::ops::Range;

use bio_types::genome::{AbstractLocus, Position};
use bio_types::strand::Strand;
use soa_derive::StructOfArray;

pub use builder::SiteMismatchesBuilder;
pub use data::{SiteData, SiteDataRef, SiteDataVec};
pub use vec::SiteMismatchesVec;

use crate::core::dna::ncounts::*;
use crate::core::dna::NucCounts;
use crate::core::dna::Nucleotide;

use super::MismatchesVec;

mod builder;
mod data;
mod vec;
