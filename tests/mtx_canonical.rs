//! Matrix Market canonicalization: duplicate coordinates are summed (the
//! Matrix Market / SciPy convention) and the resulting CSC is validated.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{ErrorCode, Reader};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_mtx_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

fn dataset(label: &str, mtx: &str) -> PathBuf {
    let d = temp_dir(label);
    write(&d.join("matrix.mtx"), mtx);
    write(
        &d.join("features.tsv"),
        "ENSG1\tA\tGene Expression\nENSG2\tB\tGene Expression\n",
    );
    write(&d.join("barcodes.tsv"), "AAAC-1\nAAAG-1\n");
    d
}

#[test]
fn duplicate_coordinates_are_summed() {
    let d = dataset(
        "dup",
        "%%MatrixMarket matrix coordinate integer general\n2 2 4\n1 1 2\n1 1 3\n2 1 1\n2 2 7\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 1, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 1.0, 7.0]);
    assert_eq!(data.metadata.stats.nnz, 3);
    assert_eq!(data.metadata.stats.total_counts, 13.0);
    assert_eq!(data.metadata.stats.min_count, 1.0);
    assert_eq!(data.metadata.stats.max_count, 7.0);
    assert_eq!(data.metadata.stats.sparsity, 0.25);
    data.matrix.validate().unwrap();
}

#[test]
fn unsorted_entries_are_canonicalized() {
    let d = dataset(
        "unsorted",
        "%%MatrixMarket matrix coordinate integer general\n2 2 3\n2 2 9\n2 1 1\n1 1 4\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 1, 1]);
    assert_eq!(data.matrix.values, vec![4.0, 1.0, 9.0]);
    data.matrix.validate().unwrap();
}

#[test]
fn duplicates_cancelling_to_zero_are_dropped() {
    let d = dataset(
        "cancel",
        "%%MatrixMarket matrix coordinate real general\n2 2 3\n1 1 2.5\n1 1 -2.5\n2 2 1\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.col_ptr, vec![0, 0, 1]);
    assert_eq!(data.matrix.row_idx, vec![1]);
    assert_eq!(data.metadata.stats.nnz, 1);
}

#[test]
fn explicit_zeros_are_not_stored() {
    let d = dataset(
        "zeros",
        "%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 0\n2 1 0\n2 2 3\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.stats.nnz, 1);
    assert_eq!(data.matrix.values, vec![3.0]);
}

#[test]
fn zero_based_index_is_rejected() {
    let d = dataset(
        "zero_idx",
        "%%MatrixMarket matrix coordinate integer general\n2 2 1\n0 1 3\n",
    );
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ValidationError);
}
