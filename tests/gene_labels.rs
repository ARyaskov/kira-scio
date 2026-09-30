//! Gene label post-processing shared by every format: Ensembl version
//! stripping (any species / feature type), duplicate detection, and the
//! opt-in scanpy-style de-duplication.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{Reader, ReaderOptions};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_labels_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

fn mtx(label: &str, features: &str, barcodes: &str) -> PathBuf {
    let d = temp_dir(label);
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 1\n3 2 2\n",
    );
    write(&d.join("features.tsv"), features);
    write(&d.join("barcodes.tsv"), barcodes);
    d
}

const VERSIONED: &str =
    "ENSDARG00000000001.3\tzf\nENSRNOG00000000002.1\trat\nENST00000380152.8\tENST00000380152.8\n";

#[test]
fn ensembl_versions_are_stripped_for_every_species_by_default() {
    let d = mtx("strip", VERSIONED, "C1\nC2\n");
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(
        data.metadata.gene_ids,
        vec![
            "ENSDARG00000000001",
            "ENSRNOG00000000002",
            "ENST00000380152"
        ]
    );
    // Symbols that are themselves Ensembl ids are stripped too.
    assert_eq!(
        data.metadata.gene_symbols,
        vec!["zf", "rat", "ENST00000380152"]
    );
}

#[test]
fn ensembl_versions_can_be_kept() {
    let d = mtx("keep", VERSIONED, "C1\nC2\n");
    let data = Reader::with_options(
        &d,
        ReaderOptions {
            strict: true,
            strip_ensembl_versions: false,
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.gene_ids[0], "ENSDARG00000000001.3");
    assert_eq!(data.metadata.gene_symbols[2], "ENST00000380152.8");
}

#[test]
fn duplicates_are_reported_for_mtx_inputs_too() {
    // Two versions of the same gene collapse to one id after stripping; the
    // duplicate is detected on the final labels.
    let d = mtx(
        "dups",
        "ENSG00000000001.1\tA\nENSG00000000001.2\tA\nENSG00000000003\tC\n",
        "C1\nC1\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(
        data.metadata.report.duplicate_gene_ids,
        vec!["ENSG00000000001"]
    );
    assert_eq!(data.metadata.report.duplicate_barcodes, vec!["C1"]);
    assert_eq!(data.metadata.gene_symbols, vec!["A", "A", "C"]);
    // Rows are kept; nothing is merged.
    assert_eq!(data.metadata.n_genes, 3);
}

#[test]
fn make_gene_labels_unique_renames_later_occurrences() {
    let d = mtx(
        "unique",
        "ENSG00000000001.1\tA\nENSG00000000001.2\tA\nENSG00000000003\tA\n",
        "C1\nC2\n",
    );
    let data = Reader::with_options(
        &d,
        ReaderOptions {
            strict: true,
            make_gene_labels_unique: true,
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(
        data.metadata.gene_ids,
        vec!["ENSG00000000001", "ENSG00000000001-1", "ENSG00000000003"]
    );
    assert_eq!(data.metadata.gene_symbols, vec!["A", "A-1", "A-2"]);
    // The report still names the original duplicates.
    assert_eq!(
        data.metadata.report.duplicate_gene_ids,
        vec!["ENSG00000000001"]
    );
}

#[test]
fn dense_duplicates_are_detected_the_same_way() {
    let d = temp_dir("dense");
    let p = d.join("m.tsv");
    write(&p, "gene\tC1\tC2\tC1\nMT-ND1\t1\t0\t2\nMT-ND1\t0\t3\t0\n");
    let data = Reader::with_options(
        &p,
        ReaderOptions {
            strict: true,
            make_gene_labels_unique: true,
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.report.duplicate_gene_ids, vec!["MT-ND1"]);
    assert_eq!(data.metadata.report.duplicate_barcodes, vec!["C1"]);
    assert_eq!(data.metadata.gene_ids, vec!["MT-ND1", "MT-ND1-1"]);
    assert_eq!(data.metadata.barcodes, vec!["C1", "C2", "C1"]);
}
