//! Matrix Market canonicalization: duplicate coordinates are summed (the
//! Matrix Market / SciPy convention) and the resulting CSC is validated.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{ErrorCode, Reader, ReaderOptions};

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

#[test]
fn entry_count_mismatch_is_an_error_in_strict_mode() {
    // Header promises 4 entries, file is truncated after 2.
    let d = dataset(
        "truncated",
        "%%MatrixMarket matrix coordinate integer general\n2 2 4\n1 1 1\n2 2 2\n",
    );
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ParseError);
    assert!(err.message.contains("declares 4"), "{}", err.message);
    assert!(err.message.contains("2 were found"), "{}", err.message);
}

#[test]
fn entry_count_mismatch_is_tolerated_in_lenient_mode() {
    let d = dataset(
        "truncated_lenient",
        "%%MatrixMarket matrix coordinate integer general\n2 2 4\n1 1 1\n2 2 2\n",
    );
    let data = Reader::with_options(
        &d,
        ReaderOptions {
            strict: false,
            force_format: None,
            feature_types: Default::default(),
            h5ad_source: Default::default(),
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.stats.nnz, 2);
}

#[test]
fn entry_count_includes_zero_and_duplicate_entries() {
    // 4 listed entries: one explicit zero, one duplicate. Count matches the
    // header even though only 2 stored values remain.
    let d = dataset(
        "count_semantics",
        "%%MatrixMarket matrix coordinate integer general\n2 2 4\n1 1 1\n1 1 2\n2 1 0\n2 2 5\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.stats.nnz, 2);
    assert_eq!(data.matrix.values, vec![3.0, 5.0]);
}

#[test]
fn header_without_entry_count_is_rejected_in_strict_mode() {
    let d = dataset(
        "no_nnz",
        "%%MatrixMarket matrix coordinate integer general\n2 2\n1 1 1\n",
    );
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ParseError);
    let data = Reader::with_options(
        &d,
        ReaderOptions {
            strict: false,
            force_format: None,
            feature_types: Default::default(),
            h5ad_source: Default::default(),
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.stats.nnz, 1);
}

/// Some exporters write Matrix Market as cells x genes. With both label
/// files present the orientation is unambiguous and the reader transposes.
#[test]
fn cells_by_genes_matrix_is_transposed_when_labels_say_so() {
    let d = temp_dir("transposed");
    // 2 rows (cells) x 3 cols (genes): cell 1 has gene1=4, gene3=1; cell 2 has gene2=7.
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 3 3\n1 1 4\n1 3 1\n2 2 7\n",
    );
    write(&d.join("features.tsv"), "ENSG1\tA\nENSG2\tB\nENSG3\tC\n");
    write(&d.join("barcodes.tsv"), "C1\nC2\n");

    let data = Reader::new(&d).read_all().unwrap();
    assert!(data.metadata.report.transposed);
    assert!(data.metadata.report.is_lossless());
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2", "ENSG3"]);
    assert_eq!(data.metadata.barcodes, vec!["C1", "C2"]);
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 2, 1]);
    assert_eq!(data.matrix.values, vec![4.0, 1.0, 7.0]);
    data.matrix.validate().unwrap();
}

#[test]
fn transposition_requires_both_label_files() {
    let d = temp_dir("transposed_one_file");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 3 1\n1 1 4\n",
    );
    write(&d.join("features.tsv"), "ENSG1\tA\nENSG2\tB\nENSG3\tC\n");
    // No barcodes.tsv: the mismatch is a plain label error in strict mode.
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::DimensionMismatch);
}

#[test]
fn correctly_oriented_matrix_is_never_transposed() {
    let d = dataset(
        "not_transposed",
        "%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 2 4\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert!(!data.metadata.report.transposed);
    assert_eq!(data.matrix.col_ptr, vec![0, 0, 1]);
    assert_eq!(data.matrix.row_idx, vec![0]);
}

/// A directory may hold several prefixed datasets; naming the matrix file
/// selects one of them, while naming the directory stays ambiguous.
#[test]
fn explicit_matrix_file_selects_its_prefix_in_a_multi_dataset_dir() {
    let d = temp_dir("multi_prefix");
    for (p, v) in [("A", "1"), ("B", "2")] {
        write(
            &d.join(format!("{p}_matrix.mtx")),
            &format!("%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 {v}\n"),
        );
        write(
            &d.join(format!("{p}_features.tsv")),
            &format!("ENSG{p}\tG{p}\n"),
        );
        write(&d.join(format!("{p}_barcodes.tsv")), &format!("C{p}\n"));
    }
    let b = Reader::new(d.join("B_matrix.mtx")).read_all().unwrap();
    assert_eq!(b.metadata.gene_ids, vec!["ENSGB"]);
    assert_eq!(b.metadata.barcodes, vec!["CB"]);
    assert_eq!(b.matrix.values, vec![2.0]);

    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ValidationError);
    assert!(err.message.contains("A, B"), "{}", err.message);
}

#[test]
fn unprefixed_matrix_file_next_to_prefixed_datasets() {
    let d = temp_dir("mixed_prefix");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 7\n",
    );
    write(&d.join("features.tsv"), "ENSG0\tG0\n");
    write(&d.join("barcodes.tsv"), "C0\n");
    write(
        &d.join("A_matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n",
    );
    write(&d.join("A_features.tsv"), "ENSGA\tGA\n");
    write(&d.join("A_barcodes.tsv"), "CA\n");
    let plain = Reader::new(d.join("matrix.mtx")).read_all().unwrap();
    assert_eq!(plain.metadata.gene_ids, vec!["ENSG0"]);
    assert_eq!(plain.matrix.values, vec![7.0]);
}

/// Matrix Market banner handling.
#[test]
fn banner_field_and_format_are_honored() {
    // `pattern`: two tokens per entry, value 1.
    let d = dataset(
        "pattern",
        "%%MatrixMarket matrix coordinate pattern general\n2 2 2\n1 1\n2 2\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.values, vec![1.0, 1.0]);
    assert!(data.metadata.stats.is_integer);

    // `real` with scientific notation and negative integer fields parse.
    let d = dataset(
        "real",
        "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1e2\n2 2 -2.5\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.values, vec![100.0, -2.5]);
    let d = dataset(
        "int_neg",
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 -3\n2 2 +4\n",
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.matrix.values, vec![-3.0, 4.0]);

    // Banner is case-insensitive.
    let d = dataset(
        "case",
        "%%matrixmarket MATRIX Coordinate Integer General\n2 2 1\n1 1 9\n",
    );
    assert_eq!(Reader::new(&d).read_all().unwrap().matrix.values, vec![9.0]);
}

#[test]
fn unsupported_banners_are_rejected() {
    for (label, banner) in [
        (
            "array",
            "%%MatrixMarket matrix array integer general\n2 2\n1\n2\n3\n4\n",
        ),
        (
            "complex",
            "%%MatrixMarket matrix coordinate complex general\n2 2 1\n1 1 1 0\n",
        ),
        (
            "symmetric",
            "%%MatrixMarket matrix coordinate integer symmetric\n2 2 1\n2 1 1\n",
        ),
        (
            "vector",
            "%%MatrixMarket vector coordinate integer general\n2 1\n1 1\n",
        ),
    ] {
        let d = dataset(label, banner);
        let err = Reader::new(&d).read_all().unwrap_err();
        assert_eq!(err.code, ErrorCode::UnsupportedFormat, "{label}: {err}");
    }
}

#[test]
fn missing_banner_is_strict_error_and_lenient_real() {
    let d = dataset("no_banner", "2 2 1\n1 1 2.5\n");
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ParseError);
    assert!(err.message.contains("banner"), "{}", err.message);
    let data = Reader::with_options(
        &d,
        ReaderOptions {
            strict: false,
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.matrix.values, vec![2.5]);
}

#[test]
fn malformed_entries_name_the_line() {
    let d = dataset(
        "malformed",
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n2 x 1\n",
    );
    let err = Reader::new(&d).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ParseError);
    assert!(err.message.contains("line 4"), "{}", err.message);
}
