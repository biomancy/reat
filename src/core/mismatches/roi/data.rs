use std::ops::Range;

use bio_types::genome::Position;
use bio_types::strand::Strand;
use soa_derive::StructOfArray;

use crate::core::dna::NucCounts;
use crate::core::mismatches::roi::NucMismatches;
use crate::core::workload::roi::*;

#[derive(Clone, Debug, StructOfArray)]
#[soa_derive(Clone, Debug)]
pub struct ROIDataRecord {
    // Original ROI record
    pub premasked: Range<Position>,
    pub postmasked: Range<Position>,
    pub subintervals: Vec<Range<Position>>,
    pub name: String,
    pub strand: Strand,
}

#[derive(Clone, Debug, StructOfArray)]
#[soa_derive(Clone, Debug)]
pub struct ROIData {
    #[nested_soa]
    pub roi: ROIDataRecord,
    // Number of unique fragments covering the ROI
    pub coverage: u32,
    // Total nucleotides in the given ROI (after masking)
    pub prednuc: NucCounts,
    // Total heterozygous loci in the ROI (after masking)
    pub heterozygous: u64,
    // Observed mismatches relative to the predicted reference
    pub mismatches: NucMismatches,
}

impl ROIDataRecordRef<'_> {
    pub fn nucmasked(&self) -> u64 {
        let mut nucin = 0;
        for piece in self.subintervals {
            nucin += piece.end - piece.start;
        }
        let total = self.premasked.end - self.premasked.start;

        total - nucin
    }
}

impl From<&'_ ROI> for ROIDataRecord {
    fn from(roi: &ROI) -> Self {
        Self {
            premasked: roi.premasked(),
            postmasked: roi.postmasked(),
            subintervals: roi.subintervals().into(),
            name: roi.name().into(),
            strand: roi.strand(),
        }
    }
}

impl From<ROIDataRecordRef<'_>> for ROIDataRecord {
    fn from(x: ROIDataRecordRef<'_>) -> Self {
        Self {
            premasked: x.premasked.to_owned(),
            postmasked: x.postmasked.to_owned(),
            subintervals: x.subintervals.to_owned(),
            name: x.name.into(),
            strand: *x.strand,
        }
    }
}

impl From<ROIDataRef<'_>> for ROIData {
    fn from(x: ROIDataRef<'_>) -> Self {
        Self {
            roi: x.roi.into(),
            coverage: *x.coverage,
            prednuc: *x.prednuc,
            heterozygous: *x.heterozygous,
            mismatches: *x.mismatches,
        }
    }
}
