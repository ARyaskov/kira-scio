//! A UTF-8 byte order mark at the start of a text file (Excel exports,
//! some Windows tools) must not leak into labels or shift columns.

use kira_scio::{DetectedFormat, Reader, detect_input_format};

mod common;
use common::{temp_dir, write};

const BOM: &str = "\u{feff}";

#[test]
fn dense_gene_major_header_with_bom() {
    let d = temp_dir("dense_gene");
    let p = d.join("m.tsv");
    write(&p, &format!("{BOM}gene\tC1\tC2\nG1\t1\t2\nG2\t0\t3\n"));
    let data = Reader::new(&p).read_all().unwrap();
    assert_eq!(data.metadata.barcodes, vec!["C1", "C2"]);
    assert_eq!(data.metadata.gene_ids, vec!["G1", "G2"]);
    assert_eq!(data.matrix.col_ptr, vec![0, 1, 3]);
    assert_eq!(data.matrix.values, vec![1.0, 2.0, 3.0]);
}

#[test]
fn dense_cell_major_header_with_bom() {
    let d = temp_dir("dense_cell");
    let p = d.join("m.csv");
    write(
        &p,
        &format!("{BOM}Cell_Index,GENE_A,GENE_B\n1,4,0\n2,0,2\n"),
    );
    let data = Reader::new(&p).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.barcodes, vec!["1", "2"]);
    assert_eq!(data.metadata.gene_symbols, vec!["GENE_A", "GENE_B"]);
}

#[test]
fn mtx_triplet_files_with_bom() {
    let d = temp_dir("mtx");
    write(
        &d.join("matrix.mtx"),
        &format!("{BOM}%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 5\n2 1 6\n"),
    );
    write(
        &d.join("features.tsv"),
        &format!("{BOM}ENSG1\tA\nENSG2\tB\n"),
    );
    write(&d.join("barcodes.tsv"), &format!("{BOM}AAAC-1\n"));
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2"]);
    assert_eq!(data.metadata.gene_symbols, vec!["A", "B"]);
    assert_eq!(data.metadata.barcodes, vec!["AAAC-1"]);
    assert_eq!(data.matrix.values, vec![5.0, 6.0]);
}

#[test]
fn sniffer_ignores_bom_when_classifying() {
    let d = temp_dir("sniff");
    let p = d.join("matrix.txt");
    write(
        &p,
        &format!("{BOM}%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 1\n"),
    );
    assert_eq!(detect_input_format(&p).unwrap(), DetectedFormat::Mtx10x);
}
