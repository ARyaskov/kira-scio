//! `IngestReport` semantics: lossy repairs are errors in strict mode and are
//! recorded in lenient mode; lossless normalizations are recorded in both.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{CountMismatch, DetectedFormat, ErrorCode, Reader, ReaderOptions};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_report_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

fn lenient(force: Option<DetectedFormat>) -> ReaderOptions {
    ReaderOptions {
        strict: false,
        force_format: force,
        feature_types: Default::default(),
        h5ad_source: Default::default(),
        ..Default::default()
    }
}

fn mtx_dataset(label: &str, mtx: &str, features: &str, barcodes: &str) -> PathBuf {
    let d = temp_dir(label);
    write(&d.join("matrix.mtx"), mtx);
    write(&d.join("features.tsv"), features);
    write(&d.join("barcodes.tsv"), barcodes);
    d
}

const FEATURES_2: &str = "ENSG1\tA\nENSG2\tB\n";
const BARCODES_2: &str = "C1\nC2\n";

#[test]
fn clean_mtx_input_yields_clean_report() {
    let d = mtx_dataset(
        "clean",
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n2 2 2\n",
        FEATURES_2,
        BARCODES_2,
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert!(data.metadata.report.is_clean());
    assert!(data.metadata.report.is_lossless());
}

#[test]
fn mtx_lenient_mode_records_every_repair() {
    // Row 5 is out of range, one NaN, one explicit zero, one duplicate, and
    // the header claims 7 entries while 6 are present. Only one feature
    // label and three barcodes are supplied for a 2x2 matrix.
    let d = mtx_dataset(
        "lenient",
        "%%MatrixMarket matrix coordinate real general\n2 2 7\n1 1 1\n5 1 1\n2 1 NaN\n2 2 0\n1 2 2\n1 2 3\n",
        "ENSG1\tA\n",
        "C1\nC2\nC3\n",
    );
    let data = Reader::with_options(&d, lenient(None)).read_all().unwrap();
    let r = &data.metadata.report;
    assert_eq!(r.dropped_out_of_range, 1);
    assert_eq!(r.dropped_non_finite, 1);
    assert_eq!(r.explicit_zeros, 1);
    assert_eq!(r.merged_duplicates, 1);
    assert_eq!(
        r.entry_count_mismatch,
        Some(CountMismatch {
            expected: 7,
            found: 6
        })
    );
    assert_eq!(
        r.relabeled_genes,
        Some(CountMismatch {
            expected: 2,
            found: 1
        })
    );
    assert_eq!(
        r.relabeled_barcodes,
        Some(CountMismatch {
            expected: 2,
            found: 3
        })
    );
    assert!(!r.is_lossless());
    assert!(!r.bom_stripped);

    // The matrix itself is what survived: (0,0)=1 and (1,0)=5.
    assert_eq!(data.matrix.col_ptr, vec![0, 1, 2]);
    assert_eq!(data.matrix.row_idx, vec![0, 0]);
    assert_eq!(data.matrix.values, vec![1.0, 5.0]);
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "gene_00000002"]);
    assert_eq!(data.metadata.barcodes, vec!["C1", "C2"]);
}

#[test]
fn mtx_strict_mode_rejects_each_lossy_repair() {
    let cases: [(&str, &str, &str, &str, ErrorCode); 4] = [
        (
            "oor",
            "%%MatrixMarket matrix coordinate integer general\n2 2 1\n5 1 1\n",
            FEATURES_2,
            BARCODES_2,
            ErrorCode::ValidationError,
        ),
        (
            "nan",
            "%%MatrixMarket matrix coordinate real general\n2 2 1\n1 1 NaN\n",
            FEATURES_2,
            BARCODES_2,
            ErrorCode::ValidationError,
        ),
        (
            "count",
            "%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 1\n",
            FEATURES_2,
            BARCODES_2,
            ErrorCode::ParseError,
        ),
        (
            "labels",
            "%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 1\n",
            "ENSG1\tA\n",
            BARCODES_2,
            ErrorCode::DimensionMismatch,
        ),
    ];
    for (label, mtx, features, barcodes, code) in cases {
        let d = mtx_dataset(label, mtx, features, barcodes);
        let err = Reader::new(&d).read_all().unwrap_err();
        assert_eq!(err.code, code, "case {label}: {err}");
    }
}

#[test]
fn mtx_lossless_normalizations_are_recorded_in_strict_mode() {
    let d = mtx_dataset(
        "lossless",
        "\u{feff}%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 1\n1 1 2\n2 2 0\n",
        FEATURES_2,
        BARCODES_2,
    );
    let data = Reader::new(&d).read_all().unwrap();
    let r = &data.metadata.report;
    assert!(r.is_lossless());
    assert!(!r.is_clean());
    assert_eq!(r.merged_duplicates, 1);
    assert_eq!(r.explicit_zeros, 1);
    assert!(r.bom_stripped);
}

#[test]
fn dense_reports_duplicate_labels_and_non_finite_values() {
    let d = temp_dir("dense");
    let p = d.join("m.tsv");
    write(
        &p,
        "gene\tC1\tC2\tC1\nMT-ND1\t1\tInf\t2\nMT-ND1\t0\t3\t0\nNDUFS1\t4\t0\t0\n",
    );
    let err = Reader::new(&p).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ValidationError);

    let data = Reader::with_options(&p, lenient(None)).read_all().unwrap();
    let r = &data.metadata.report;
    assert_eq!(r.dropped_non_finite, 1);
    assert_eq!(r.duplicate_gene_ids, vec!["MT-ND1"]);
    assert_eq!(r.duplicate_barcodes, vec!["C1"]);
    assert_eq!(r.explicit_zeros, 0, "dense inputs do not count zeros");
    assert_eq!(r.merged_duplicates, 0);
    assert!(!r.is_lossless());
    // Duplicate labels keep their own rows/columns.
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(data.metadata.n_cells, 3);
}

#[test]
fn dense_clean_input_is_clean_and_report_survives_bd_shim() {
    let d = temp_dir("bd");
    let p = d.join("S_RSEC_MolsPerCell.csv");
    write(&p, "\u{feff}Cell_Index,G1,G2\n1,1,0\n2,0,2\n");
    let md = Reader::new(&p).read_metadata().unwrap();
    assert_eq!(md.format, "bd_rhapsody_wta");
    assert!(md.report.is_lossless());
    assert!(md.report.bom_stripped);
}
