pub use buffers::{CountsBuffer, CountsBufferContent, LocusCounts, StrandedCountsBuffer, UnstrandedCountsBuffer};
pub use nuccounter::{BaseNucCounter, NucCounter, NucCounterContent};

#[cfg(test)]
pub use nuccounter::MockNucCounter;

mod buffers;
mod nuccounter;
