//! BD Rhapsody Sequence Analysis Pipeline output: `*_MolsPerCell.csv` tables
//! are cell-major (one row per cell, `Cell_Index` first column) and start
//! with `#` comment lines.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{DetectedFormat, Reader, detect_input_format};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_bd_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

const MOLS_PER_CELL: &str = "\
####################\n\
## BD Rhapsody Sequence Analysis Pipeline Version 2.0\n\
####################\n\
Cell_Index,GENE_A,GENE_B,GENE_C\n\
1,5,0,1\n\
2,0,3,0\n";

#[test]
fn mols_per_cell_csv_is_cell_major() {
    let d = temp_dir("orientation");
    let p = d.join("Sample_RSEC_MolsPerCell.csv");
    write(&p, MOLS_PER_CELL);

    assert_eq!(
        detect_input_format(&p).unwrap(),
        DetectedFormat::BdRhapsodyWta
    );
    let data = Reader::new(&p).read_all().unwrap();
    assert_eq!(data.metadata.format, "bd_rhapsody_wta");
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(data.metadata.barcodes, vec!["1", "2"]);
    assert_eq!(
        data.metadata.gene_symbols,
        vec!["GENE_A", "GENE_B", "GENE_C"]
    );
    // Cell 1: GENE_A=5, GENE_C=1; cell 2: GENE_B=3.
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 2, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 1.0, 3.0]);
    assert_eq!(data.metadata.stats.nnz, 3);
    assert_eq!(data.metadata.stats.total_counts, 9.0);
}

#[test]
fn pipeline_directory_prefers_dbec_over_rsec() {
    let d = temp_dir("dir_dbec");
    write(&d.join("Sample_RSEC_MolsPerCell.csv"), MOLS_PER_CELL);
    write(
        &d.join("Sample_DBEC_MolsPerCell.csv"),
        "Cell_Index,GENE_A,GENE_B\n1,4,0\n2,0,2\n3,1,1\n",
    );
    write(
        &d.join("Sample_Bioproduct_Stats.csv"),
        "Bioproduct,Reads\nGENE_A,10\n",
    );

    assert_eq!(
        detect_input_format(&d).unwrap(),
        DetectedFormat::BdRhapsodyWta
    );
    let resolved = kira_scio::resolve_bd_input_path(&d).unwrap();
    assert_eq!(
        resolved.file_name().unwrap().to_str().unwrap(),
        "Sample_DBEC_MolsPerCell.csv"
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 3);
    assert_eq!(data.metadata.n_genes, 2);
}

#[test]
fn pipeline_directory_falls_back_to_rsec() {
    let d = temp_dir("dir_rsec");
    write(&d.join("Sample_RSEC_MolsPerCell.csv"), MOLS_PER_CELL);
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 3);
}

#[test]
fn legacy_raw_counts_still_wins_over_mols_per_cell() {
    let d = temp_dir("dir_legacy");
    write(&d.join("Sample_RSEC_MolsPerCell.csv"), MOLS_PER_CELL);
    write(&d.join("raw_counts.tsv"), "gene\tC1\nG1\t7\n");
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 1);
    assert_eq!(data.metadata.barcodes, vec!["C1"]);
}

#[test]
fn mex_export_directory_reads_as_mtx() {
    // `<sample>_RSEC_MolsPerCell_MEX/` uses the 10x MEX triplet.
    let d = temp_dir("mex").join("Sample_RSEC_MolsPerCell_MEX");
    fs::create_dir_all(&d).unwrap();
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 3\n2 2 4\n",
    );
    write(
        &d.join("features.tsv"),
        "GENE_A\tGENE_A\tGene Expression\nGENE_B\tGENE_B\tGene Expression\n",
    );
    write(&d.join("barcodes.tsv"), "1\n2\n");
    assert_eq!(detect_input_format(&d).unwrap(), DetectedFormat::Mtx10x);
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.gene_symbols, vec!["GENE_A", "GENE_B"]);
}
