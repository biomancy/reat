pub use buffers::{CountsBuffer, CountsBufferContent, NucCounts, StrandedCountsBuffer, UnstrandedCountsBuffer};
pub use nuccounter::{BaseNucCounter, NucCounter, NucCounterContent};

#[cfg(test)]
pub use nuccounter::MockNucCounter;

mod buffers;
mod nuccounter;
