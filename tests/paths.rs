use std::path::PathBuf;

use static_init::dynamic;

#[dynamic]
pub static RESOURCES: PathBuf = PathBuf::from(file!()).parent().unwrap().join("resources");

#[allow(non_snake_case)]
pub mod GRCh38 {
    use super::*;

    #[dynamic]
    pub static FOLDER: PathBuf = RESOURCES.join("GRCh38");

    #[dynamic]
    pub static GFF3: String = FOLDER.join("Homo_sapiens.GRCh38.104.gff3.gz").to_str().unwrap().to_string();

    #[dynamic]
    pub static FASTA: String =
        FOLDER.join("Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz").to_str().unwrap().to_string();

    #[dynamic]
    pub static ALU: String = FOLDER.join("Alu.bed.gz").to_str().unwrap().to_string();
}

pub mod bam {
    use super::*;

    #[dynamic]
    pub static FOLDER: PathBuf = RESOURCES.join("bam");

    #[dynamic]
    pub static EXAMPLE: String = FOLDER.join("SRX6966474.bam").to_str().unwrap().to_string();
}

pub mod expected {
    use super::*;

    #[dynamic]
    pub static FOLDER: PathBuf = RESOURCES.join("expected");

    #[dynamic]
    pub static LOCI: PathBuf = FOLDER.join("loci");

    #[dynamic]
    pub static ROIS: PathBuf = FOLDER.join("rois");
}
