use std::collections::HashMap;
use std::path::Path;

use bio_types::genome::Position;
use itertools::Itertools;
use rust_htslib::bcf::header::Id;
use rust_htslib::bcf::record::GenotypeAllele;
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

        let mut genotype = record
            .genotypes()
            .unwrap()
            .get(0)
            .iter()
            .filter_map(|g| match g {
                GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => None,
                GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => Some(*i as usize),
            })
            .collect_vec();
        // Unknown genotype + known -> set to known only
        if genotype.len() == 1 {
            genotype.push(genotype[0])
        }
        // Skip reference or strange genotypes
        if genotype.len() == 0 || (genotype[0] == 0 && genotype[1] == 0) || genotype.len() > 2 {
            continue;
        }

        // Fetch alleles
        let (first, second) = (record.alleles()[genotype[0]], record.alleles()[genotype[1]]);

        // Skip indels or large variants
        if first.len() != 1 || second.len() != 1 {
            continue;
        }
        let (first, second): (ReqNucleotide, ReqNucleotide) = match (first[0].try_into(), second[0].try_into()) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };

        let rid = record.rid().unwrap() as usize;
        match first == second {
            true => homozygous[rid].push((record.pos() as Position, first)),
            false => heterozygous[rid].push((record.pos() as Position, first, second)),
        }
    }
    SimplisticSNV { rid2ref, ref2rid, homozygous, heterozygous }
}
