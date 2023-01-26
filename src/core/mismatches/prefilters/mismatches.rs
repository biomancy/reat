use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{ROIData, ROINucCounts};
use crate::core::mismatches::site::SiteData;
use crate::core::refpred::PredNucleotide;

use super::MismatchesPreFilter;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct ByMismatches {
    minfreq: f32,
    // Precasted values to save on convertions
    minmismatches_f32: f32,
    mincov_f32: f32,
    minmismatches_u32: u32,
    mincov_u32: u32,
}

impl ByMismatches {
    pub fn new(minmismatches: u32, minfreq: f32, mincov: u32) -> Self {
        Self {
            minfreq,
            minmismatches_f32: minmismatches as f32,
            mincov_f32: mincov as f32,
            minmismatches_u32: minmismatches,
            mincov_u32: mincov,
        }
    }

    #[inline]
    pub fn enough_mismatches_per_roi(&self, x: &ROINucCounts) -> bool {
        let (cov, mismatch) = (x.coverage(), x.mismatches());
        cov >= self.mincov_f32 && mismatch >= self.minmismatches_f32 && mismatch / cov >= self.minfreq
    }

    #[inline]
    pub fn enough_mismatches_per_site(&self, reference: Nucleotide, sequenced: &NucCounts) -> bool {
        let cov = sequenced.coverage();
        let mismatch = sequenced.mismatches(reference);
        cov >= self.mincov_u32 && mismatch >= self.minmismatches_u32 && mismatch as f32 / cov as f32 >= self.minfreq
    }

    #[inline]
    pub fn mincov(&self) -> u32 {
        self.mincov_u32
    }

    #[inline]
    pub fn minfreq(&self) -> f32 {
        self.minfreq
    }

    #[inline]
    pub fn minmismatches(&self) -> u32 {
        self.minmismatches_u32
    }
}

impl MismatchesPreFilter<ROIData> for ByMismatches {
    #[inline]
    fn is_ok(&self, preview: &ROIData) -> bool {
        self.enough_mismatches_per_roi(&preview.mismatches)
    }
}

impl MismatchesPreFilter<SiteData> for ByMismatches {
    #[inline]
    fn is_ok(&self, preview: &SiteData) -> bool {
        match preview.prednuc {
            PredNucleotide::Homozygous(nuc) => self.enough_mismatches_per_site(nuc, &preview.sequenced),
            PredNucleotide::Heterozygous((n1, n2)) => {
                self.enough_mismatches_per_site(n1, &preview.sequenced)
                    || self.enough_mismatches_per_site(n2, &preview.sequenced)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::dna::{NucCounts, Nucleotide};
    use crate::core::mismatches::roi::ROINucCounts;

    use super::*;

    #[test]
    fn ok_roi() {
        let mut dummy: ROINucCounts = Default::default();

        dummy.A.C = 1_f32;
        dummy.A.A = 4_f32;

        dummy.C.G = 1_f32;
        dummy.C.C = 2_f32;

        dummy.G.T = 1_f32;
        dummy.G.G = 5_f32;

        dummy.T.A = 10_f32;
        dummy.T.T = 3_f32;
        // dummy coverage = 27, mismatches = 13, freq = 0.48148

        for (expected, minmismatches, minfreq, mincov) in [
            (false, 14, 0f32, 0),
            (true, 13, 0f32, 0),
            (true, 12, 0f32, 0),
            (true, 13, 0.48f32, 0),
            (false, 13, 0.5f32, 0),
            (true, 13, 0.48f32, 11),
            (true, 13, 0.48f32, 27),
            (false, 13, 0.48f32, 30),
            (false, 1, 0.48f32, 42),
        ] {
            let filter = ByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(filter.enough_mismatches_per_roi(&dummy), expected, "{} {} {}", minmismatches, minfreq, mincov);
        }
    }

    #[test]
    fn ok_site() {
        let mut reference = Nucleotide::A;
        let sequenced = NucCounts { A: 1, C: 2, G: 3, T: 4 };

        for (expected, minmismatches, minfreq, mincov) in [
            (false, 10, 0f32, 0),
            (true, 9, 0f32, 5),
            (true, 8, 0f32, 8),
            (true, 9, 0.85f32, 9),
            (false, 9, 0.95f32, 10),
            (true, 9, 0.85f32, 10),
            (false, 9, 0.85f32, 11),
        ] {
            let filter = ByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(filter.enough_mismatches_per_site(reference, &sequenced), expected);
        }

        reference = Nucleotide::Unknown;
        for (expected, minmismatches, minfreq, mincov) in [
            (true, 10, 0f32, 0),
            (true, 9, 0f32, 0),
            (false, 11, 0f32, 0),
            (true, 10, 1f32, 0),
            (false, 11, 1f32, 0),
            (true, 10, 1f32, 10),
            (false, 10, 1f32, 11),
        ] {
            let filter = ByMismatches::new(minmismatches, minfreq, mincov);
            assert_eq!(filter.enough_mismatches_per_site(reference, &sequenced), expected);
        }
    }
}
