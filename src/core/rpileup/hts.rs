use std::path::PathBuf;

use bio_types::genome::AbstractInterval;
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read, Record};

use crate::core::rpileup::{ReadsCollider, ReadsCollidingEngine};

pub struct HTSPileupEngine<Collider> {
    collider: Collider,
    htsreaders: Vec<IndexedReader>,
    htsfiles: Vec<PathBuf>,
    success: bool,
}

impl<Collider: for<'a> ReadsCollider<'a, Record>> HTSPileupEngine<Collider> {
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

        Self { collider, htsreaders, htsfiles, success: false }
    }
}

impl<Collider: for<'a> ReadsCollider<'a, Record>> ReadsCollidingEngine<Record, Collider> for HTSPileupEngine<Collider> {
    fn run(&mut self, cwork: <Collider as ReadsCollider<'_, Record>>::Workload) {
        let toread = self
            .htsreaders
            .iter_mut()
            .filter(|x| x.header().target_names().contains(&cwork.contig().as_bytes()))
            .map(|x| {
                x.fetch((cwork.contig(), cwork.range().start, cwork.range().end)).unwrap_or_else(|_| {
                    panic!(
                        "Failed to fetch reads for {}:{}-{} (HTS file corrupted?)",
                        cwork.contig(),
                        cwork.range().start,
                        cwork.range().end
                    )
                });

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

    fn result(&self) -> Result<<Collider as ReadsCollider<'_, Record>>::ColliderResult, ()> {
        if self.success {
            Ok(self.collider.result())
        } else {
            Err(())
        }
    }
}

impl<Collider: for<'a> ReadsCollider<'a, Record> + Clone> Clone for HTSPileupEngine<Collider> {
    fn clone(&self) -> Self {
        Self::new(self.htsfiles.clone(), self.collider.clone())
    }
}
