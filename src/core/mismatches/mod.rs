use bio_types::strand::Strand;

pub mod interval;
pub mod roi;
pub mod site;

pub trait IntermediateMismatches {
    type Finalized;

    fn finalize(self) -> Self::Finalized;
    fn set_strand(&mut self, strand: Strand);
}
