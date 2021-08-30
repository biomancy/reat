use std::path::Path;
use std::path::PathBuf;

const RESOURCES: PathBuf = PathBuf::from(file!()).parent().unwrap().join("resources");

#[allow(non_snake_case)]
pub mod GRCh38 {
    use super::*;

    const FOLDER: PathBuf = RESOURCES.join("GRCh38");

    pub fn gtf() -> &'static Path {
        &FOLDER.join("Homo_sapiens.GRCh38.104.gtf.gz")
    }
}
