use std::path::PathBuf;

use bio_types::genome::{AbstractInterval, Interval};
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read, Record};

use crate::core::read::AlignedRead;
use crate::core::rpileup::{ReadsCollider, ReadsPileupEngine};
use std::marker::PhantomData;

pub struct HTSPileupEngine<'a, Collider: ReadsCollider<'a, Record>> {
    collider: Collider,
    htsreaders: Vec<IndexedReader>,
    htsfiles: Vec<PathBuf>,
    success: bool,
    marker: PhantomData<&'a ()>,
}

impl<'a, Collider: ReadsCollider<'a, Record>> HTSPileupEngine<'a, Collider> {
    pub fn new(htsfiles: Vec<PathBuf>, collider: Collider) -> Self {
        let htsreaders: Vec<IndexedReader> = htsfiles
            .iter()
            .map(|hts| {
                IndexedReader::from_path(&hts).unwrap_or_else(|_| {
                    panic!(
                        "Failed to open file {}\n\
                        Possible reasons: BAM file was not indexed (samtools index); you don't have read permissions",
                        hts.display()
                    )
                })
            })
            .collect();

        Self { collider, htsreaders, htsfiles, success: false, marker: Default::default() }
    }
}

impl<'a, Collider: ReadsCollider<'a, Record>> ReadsPileupEngine<'a, Record, Collider>
    for HTSPileupEngine<'a, Collider>
{
    fn run(&mut self, interval: Interval, cwork: Collider::Workload) {
        let toread = self
            .htsreaders
            .iter_mut()
            .filter(|x| x.header().target_names().contains(&interval.contig().as_bytes()))
            .map(|x| {
                x.fetch((interval.contig(), interval.range().start, interval.range().end))
                    .unwrap_or_else(|_| panic!("Failed to fetch filters for {:?} (HTS file corrupted?)", interval));

                let mut record = Record::new();
                let r = x.read(&mut record);
                (x, record, r)
            })
            .filter(|(_, _, r)| r.as_ref().map(|x| x.is_ok()).unwrap_or(false))
            .map(|(x, record, r)| {
                debug_assert!(r.is_some() && r.unwrap().is_ok());
                (x, record)
            })
            .collect_vec();

        // Nothing to do
        if toread.is_empty() {
            self.success = false;
            return;
        }

        // Something to do, trigger the reset
        self.collider.reset(cwork);

        for (reader, mut record) in toread.into_iter() {
            self.collider.collide(&record);
            while let Some(Ok(())) = reader.read(&mut record) {
                self.collider.collide(&record);
            }
        }
        self.collider.finalize();
        self.success = true;
    }

    fn result(&'a self) -> Result<Collider::ColliderResult, ()> {
        if self.success {
            Ok(self.collider.result())
        } else {
            Err(())
        }
    }
}

impl<'a, Collider: ReadsCollider<'a, Record> + Clone> Clone for HTSPileupEngine<'a, Collider> {
    fn clone(&self) -> Self {
        Self::new(self.htsfiles.clone(), self.collider.clone())
    }
}
