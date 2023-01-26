use std::collections::HashMap;
use std::ops::Range;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::genome::Position;
use bio_types::strand::Strand;

use crate::core::dna::NucCounts;
use crate::core::dna::ReqNucleotide;
use crate::core::io::fasta::FastaReader;
use crate::core::io::vcf;
use crate::core::refpred::PredNucleotide::{Heterozygous, Homozygous};
use crate::core::refpred::{PredNucleotide, RefEngineResult};

use super::RefEngine;

#[derive(Clone)]
pub struct VCFCorrectedReference {
    ref2rid: HashMap<String, u32>,
    homozygotes: AnnotMap<u32, ReqNucleotide>,
    heterozygotes: AnnotMap<u32, (ReqNucleotide, ReqNucleotide)>,
    cache: Vec<PredNucleotide>,
    reader: Box<dyn FastaReader>,
}

impl VCFCorrectedReference {
    pub fn new(variants: vcf::SimplisticSNV, reader: Box<dyn FastaReader>) -> Self {
        debug_assert_eq!(variants.homozygous.len(), variants.rid2ref.len());
        debug_assert_eq!(variants.heterozygous.len(), variants.rid2ref.len());
        debug_assert_eq!(variants.ref2rid.len(), variants.rid2ref.len());
        // Parse homozygotes
        let mut homozygotes = AnnotMap::new();
        for (rid, loci) in variants.homozygous.into_iter().enumerate() {
            for (locus, nuc) in loci {
                let record = Contig::new(rid as u32, locus as isize, 1, Strand::Unknown);
                homozygotes.insert_at(nuc, &record);
            }
        }
        // Parse heterozygotes
        let mut heterozygotes = AnnotMap::new();
        for (rid, loci) in variants.heterozygous.into_iter().enumerate() {
            for (locus, nuc1, nuc2) in loci {
                let record = Contig::new(rid as u32, locus as isize, 1, Strand::Unknown);
                heterozygotes.insert_at((nuc1, nuc2), &record);
            }
        }
        Self { ref2rid: variants.ref2rid, homozygotes, heterozygotes, cache: vec![], reader }
    }
}

impl RefEngine for VCFCorrectedReference {
    fn run(&mut self, contig: &str, range: Range<Position>, _: &[NucCounts]) {
        let length = (range.end - range.start) as usize;
        self.cache.clear();
        self.cache.reserve(length);

        self.reader.fetch(contig, range.clone());
        let reference = self.reader.result();

        for refn in reference.iter() {
            self.cache.push(Homozygous(*refn));
        }

        let rid = *self.ref2rid.get(contig).expect("Failed to fetch VCF data variants for a contig.");
        let query = Contig::new(rid, range.start as isize, length, Strand::Unknown);
        // Fill in homozygous positions
        for hit in self.homozygotes.find(&query) {
            let loc = hit.interval().start as Position;
            debug_assert!(range.contains(&loc));
            self.cache[(loc - range.start) as usize] = Homozygous((*hit.data()).into())
        }
        // Fill in heterozygous positions
        for hit in self.heterozygotes.find(&query) {
            let loc = hit.interval().start as Position;
            debug_assert!(range.contains(&loc));
            let (f, s) = hit.data();
            self.cache[(loc - range.start) as usize] = Heterozygous(((*f).into(), (*s).into()))
        }
    }

    fn results(&self) -> RefEngineResult<'_> {
        RefEngineResult { predicted: &self.cache, reference: self.reader.result() }
    }
}
