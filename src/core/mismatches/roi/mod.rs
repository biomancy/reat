use std::ops::Range;

use bio_types::genome::{AbstractInterval, Position};
use bio_types::strand::Strand;
use soa_derive::StructOfArray;

pub use builder::ROIMismatchesBuilder;
pub use data::{ROIData, ROIDataRecord, ROIDataRecordRef, ROIDataRecordVec, ROIDataRef, ROIDataVec};
pub use msummary::NucMismatches;
use msummary::*;
pub use vec::ROIMismatchesVec;

use crate::core::dna::ncounts::*;
use crate::core::dna::NucCounts;
use crate::core::mismatches::MismatchesVec;
use crate::core::workload::roi::*;

mod builder;
mod data;
mod msummary;
mod vec;
