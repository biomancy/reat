use std::collections::HashMap;
use std::path::Path;

use bio_types::genome::Position;
use rust_htslib::bcf::header::Id;
use rust_htslib::bcf::{Read, Reader};

use crate::core::dna::ReqNucleotide;

#[derive(Default)]
pub struct SimplisticSNV {
    pub rid2ref: Vec<String>,
    pub ref2rid: HashMap<String, u32>,
    pub homozygous: Vec<Vec<(Position, ReqNucleotide)>>,
    pub heterozygous: Vec<Vec<(Position, ReqNucleotide, ReqNucleotide)>>,
}

pub fn parse(vcf: impl AsRef<Path>) -> SimplisticSNV {
    let mut reader = Reader::from_path(vcf).expect("Error opening file.");

    // Header rid2ref / ref2rid
    let contigs = reader.header().contig_count() as usize;
    let mut rid2ref = Vec::with_capacity(contigs);
    let mut ref2rid = HashMap::with_capacity(contigs);
    for contig in 0..contigs {
        let contig = contig as u32;
        let name = reader.header().rid2name(contig).expect("Failed to parse VCF contig names");
        rid2ref.push(String::from_utf8_lossy(name).into_owned());
        ref2rid.insert(String::from_utf8_lossy(name).into_owned(), contig);
    }
    // Parsed variants
    let (mut homozygous, mut heterozygous) = (vec![Vec::new(); contigs], vec![Vec::new(); contigs]);

    let mut record = reader.empty_record();
    while let Some(Ok(())) = reader.read(&mut record) {
        assert!(record.pos() >= 0, "There are negative positions in the VCF file");
        assert_eq!(record.sample_count(), 1, "Each vcf must contain exactly 1 sample");

        // Consider only records that passed all filters
        if !record.has_filter(&Id(0)) {
            continue;
        }

        let alleles = record.alleles();
        // Skip complex mismatches
        if alleles.len() != 2 || alleles[0].len() != 1 || alleles[1].len() != 1 {
            continue;
        }

        // Consider only well-defined mismatches
        let (refn, altn): (ReqNucleotide, ReqNucleotide) = match (alleles[0][0].try_into(), alleles[1][0].try_into()) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };

        // First/Second locus
        let genotype = record.genotypes().expect("Failed to read genotypes information").get(0);
        let (floc, sloc) = (genotype[0].index(), genotype[1].index());
        if floc.is_none() || sloc.is_none() {
            continue;
        }
        let (floc, sloc) = (floc.unwrap(), sloc.unwrap());
        match (floc, sloc) {
            // Reference genotype => do nothing
            (0, 0) => continue,
            // Heterozygous
            (0, 1) | (1, 0) => {
                let rid = record.rid().unwrap() as usize;
                heterozygous[rid].push((record.pos() as Position, refn, altn));
            }
            // Alt genotype
            (1, 1) => {
                let rid = record.rid().unwrap() as usize;
                homozygous[rid].push((record.pos() as Position, altn));
            }
            _ => unreachable!("Only two alleles available, but the sample genotypes contains values > 1"),
        }
    }
    SimplisticSNV { rid2ref, ref2rid, homozygous, heterozygous }
}
