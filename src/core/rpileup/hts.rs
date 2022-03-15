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
            .filter_map(|reader| {
                // No such contig in the BAM file
                if !reader.header().target_names().contains(&cwork.contig().as_bytes()) {
                    return None;
                };

                reader.fetch((cwork.contig(), cwork.range().start, cwork.range().end)).unwrap_or_else(|_| {
                    panic!(
                        "Failed to fetch reads for {}:{}-{} (HTS file corrupted?)",
                        cwork.contig(),
                        cwork.range().start,
                        cwork.range().end
                    )
                });

                let mut record = Record::new();
                let isavailable = reader.read(&mut record);
                match isavailable {
                    Some(Ok(())) => Some((reader, record)),
                    _ => None,
                }
            })
            .collect_vec();

        // Nothing to do
        if toread.is_empty() {
            self.success = false;
            return;
        }

        // Something to do, trigger the reset -> collide -> finalize
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

    fn result(&self) -> Option<<Collider as ReadsCollider<'_, Record>>::ColliderResult> {
        match self.success {
            true => Some(self.collider.result()),
            false => None,
        }
    }
}

impl<Collider: for<'a> ReadsCollider<'a, Record> + Clone> Clone for HTSPileupEngine<Collider> {
    fn clone(&self) -> Self {
        Self::new(self.htsfiles.clone(), self.collider.clone())
    }
}
