pub mod records;
pub mod summary;

// pub struct ByMismatches {
//     pub freq: f32,
//     pub number: u32
// }
//
// pub struct ByRegions {
//
// }
//
// pub struct Filters {
//     by_quality: BySequencingQuality,
//     by_mismatches: Option<ByMismatches>,
//     by_regions: Option<ByRegions>,
// }
//
// impl Filters {
//     pub fn by_quality(&self) -> &BySequencingQuality {
//         &self.by_quality
//     }
//
//     pub fn by_mismatches(&self) -> &Option<ByMismatches> {
//         &self.by_mismatches
//     }
//
//     pub fn by_regions(&self) -> &Option<ByRegions> {
//         &self.by_regions
//     }
// }
//
//
// pub struct FiltersBuilder {
//     by_quality: BySequencingQuality,
//     by_mismatches: Option<ByMismatches>,
//     by_regions: Option<ByRegions>
// }
//
// impl FiltersBuilder {
//     pub fn new() -> FiltersBuilder {
//         let by_quality =
//         FiltersBuilder { by_quality, by_mismatches: None, by_regions: None }
//     }
//
//     pub fn by_quality(mut self, by_quality: BySequencingQuality) -> FiltersBuilder {
//         self.by_quality = by_quality;
//         self
//     }
//
//     pub fn by_mismatches(mut self, by_mismatches: ByMismatches) -> FiltersBuilder {
//         self.by_mismatches = Some(by_mismatches);
//         self
//     }
//
//     pub fn by_regions(mut self, by_regions: ByRegions) -> FiltersBuilder {
//         self.by_regions = Some(by_regions);
//         self
//     }
//
//     pub fn get(self) -> Filters {
//         Filters {by_quality: self.by_quality, by_mismatches: self.by_mismatches, by_regions: self.by_regions}
//     }
// }
//
//
