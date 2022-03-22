use derive_getters::Getters;
use derive_more::Constructor;

use crate::core::dna::{NucCounts, Nucleotide};
use crate::core::mismatches::roi::{NucMismatches, ROIData};
use crate::core::mismatches::site::SiteData;
use crate::core::refpred::PredNucleotide;

use super::MismatchesPreFilter;

#[derive(Constructor, Getters, Debug, PartialEq, Copy, Clone)]
pub struct ByMismatches {
    minmismatches: u32,
    minfreq: f32,
    mincov: u32,
}

impl ByMismatches {
    #[inline]
    pub fn enough_mismatches_per_roi(&self, x: &NucMismatches) -> bool {
        let (cov, mismatch) = (x.coverage(), x.mismatches());
        cov >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
    }

    #[inline]
    pub fn enough_mismatches_per_site(&self, reference: Nucleotide, sequenced: &NucCounts) -> bool {
        let cov = sequenced.coverage();
        let mismatch = sequenced.mismatches(reference);
        cov >= self.mincov && mismatch >= self.minmismatches && mismatch as f32 / cov as f32 >= self.minfreq
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
    use crate::core::mismatches::roi::NucMismatches;

    use super::*;

    #[test]
    fn ok_roi() {
        let mut dummy: NucMismatches = Default::default();

        dummy.A.C = 1;
        dummy.A.A = 4;

        dummy.C.G = 1;
        dummy.C.C = 2;

        dummy.G.T = 1;
        dummy.G.G = 5;

        dummy.T.A = 10;
        dummy.T.T = 3;
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
            let filter = ByMismatches { minmismatches, minfreq, mincov };
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
            let filter = ByMismatches { minmismatches, minfreq, mincov };
            assert_eq!(filter.enough_mismatches_per_site(reference, &sequenced), expected);
        }
    }
}
