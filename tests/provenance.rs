//! `InputMetadata::provenance`: resolved source file, format dialect,
//! dataset prefix and raw/filtered kind inferred from naming.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{MatrixKind, Reader};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_prov_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

const MTX: &str = "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n";

#[test]
fn mtx_v3_prefixed_raw_directory() {
    let d = temp_dir("mtxv3").join("raw_feature_bc_matrix");
    fs::create_dir_all(&d).unwrap();
    write(&d.join("S1_matrix.mtx"), MTX);
    write(&d.join("S1_features.tsv"), "ENSG1\tA\tGene Expression\n");
    write(&d.join("S1_barcodes.tsv"), "C1\n");
    let md = Reader::new(&d).read_metadata().unwrap();
    let p = &md.provenance;
    assert_eq!(p.dialect, "mtx-v3");
    assert_eq!(p.dataset_prefix.as_deref(), Some("S1"));
    assert_eq!(p.matrix_kind, MatrixKind::Raw);
    assert_eq!(p.source_path, d.join("S1_matrix.mtx"));
    assert_eq!(p.matrix_source, None);
}

#[test]
fn mtx_v2_filtered_directory() {
    let d = temp_dir("mtxv2")
        .join("filtered_gene_bc_matrices")
        .join("GRCh38");
    fs::create_dir_all(&d).unwrap();
    write(&d.join("matrix.mtx"), MTX);
    write(&d.join("genes.tsv"), "ENSG1\tA\n");
    write(&d.join("barcodes.tsv"), "C1\n");
    let md = Reader::new(&d).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "mtx-v2");
    assert_eq!(md.provenance.dataset_prefix, None);
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Filtered);
}

#[test]
fn dense_dialects_and_unknown_kind() {
    let d = temp_dir("dense");
    let g = d.join("counts.tsv");
    write(&g, "gene\tC1\nG1\t1\n");
    let md = Reader::new(&g).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "dense-gene-major");
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Unknown);
    assert_eq!(md.provenance.source_path, g);

    let c = d.join("counts.csv");
    write(&c, "cell,G1\nC1,1\n");
    let md = Reader::new(&c).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "dense-cell-major");
}

#[test]
fn bd_rhapsody_dialects_and_kinds() {
    let d = temp_dir("bd");
    let dbec = d.join("S_DBEC_MolsPerCell.csv");
    write(&dbec, "Cell_Index,G1\n1,2\n");
    let md = Reader::new(&dbec).read_metadata().unwrap();
    assert_eq!(md.format, "bd_rhapsody_wta");
    assert_eq!(md.provenance.dialect, "bd-molspercell-dbec");
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Filtered);
    assert_eq!(md.provenance.source_path, dbec);

    let unf = d.join("S_RSEC_MolsPerCell_Unfiltered.csv");
    write(&unf, "Cell_Index,G1\n1,2\n");
    let md = Reader::new(&unf).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "bd-molspercell-rsec");
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Raw);

    let legacy = d.join("raw_counts.tsv");
    write(&legacy, "gene\tC1\nG1\t1\n");
    let md = Reader::new(&d).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "bd-raw-counts");
    assert_eq!(md.provenance.source_path, legacy);
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Unknown);
}
