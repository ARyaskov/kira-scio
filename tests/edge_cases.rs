use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{DetectedFormat, Reader, ReaderOptions, detect_input_format};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_edge_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &PathBuf, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

fn temp_file(label: &str, ext: &str, content: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("kira_scio_edge_{label}_{ts}.{ext}"));
    write(&path, content);
    path
}

/// `gene_ids` and `gene_symbols` should expose distinct columns.
#[test]
fn s1_features_tsv_separates_id_and_symbol() {
    let d = temp_dir("s1");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate real general\n2 1 1\n1 1 1.0\n",
    );
    write(&d.join("features.tsv"), "ENSG1\tBRCA1\nENSG2\tTP53\n");
    write(&d.join("barcodes.tsv"), "C1\n");

    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2"]);
    assert_eq!(data.metadata.gene_symbols, vec!["BRCA1", "TP53"]);
}

/// duplicate gene symbols are kept in insertion order rather than
/// triggering a hard failure in strict mode.
#[test]
fn s2_duplicate_genes_are_kept_in_strict() {
    let p = temp_file(
        "s2",
        "tsv",
        "#meta\ncellA\tcellB\nMT-ND1\t1\t2\nMT-ND1\t3\t4\nNDUFS1\t5\t6\n",
    );
    let data = Reader::with_options(
        &p,
        ReaderOptions {
            strict: true,
            force_format: Some(DetectedFormat::BdRhapsodyWta),
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.gene_symbols, vec!["MT-ND1", "MT-ND1", "NDUFS1"]);
}

/// content sniffing should classify a `.txt` file with a leading comment
/// and float values as BD Rhapsody.
#[test]
fn s3_content_sniff_promotes_bd_for_floats() {
    let p = temp_file("s3a", "txt", "#meta\ncellA\tcellB\nGENE1\t1.0\t2.5\n");
    assert_eq!(detect_input_format(&p).unwrap(), DetectedFormat::BdRhapsodyWta);
}

/// a plain integer-only `.txt` should default to DenseTsvCsv, **not** BD.
#[test]
fn s3_content_sniff_keeps_integers_as_dense() {
    let p = temp_file("s3b", "txt", "cellA\tcellB\nGENE1\t1\t2\n");
    assert_eq!(detect_input_format(&p).unwrap(), DetectedFormat::DenseTsvCsv);
}

/// file whose name contains `_raw_counts.tsv` (anywhere) should resolve
/// to BD Rhapsody, even if the actual extension is `.txt`.
#[test]
fn s3_filename_substring_matches_bd() {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let p = std::env::temp_dir().join(format!("sample_raw_counts.tsv_{ts}.txt"));
    write(&p, "gene\tcellA\tcellB\nGENE1\t1\t2\n");
    assert_eq!(detect_input_format(&p).unwrap(), DetectedFormat::BdRhapsodyWta);
}

/// dense parser must treat a leading empty header cell (`\tC1\tC2`) as a
/// row-major matrix.
#[test]
fn s4_empty_first_header_cell() {
    let p = temp_file("s4", "tsv", "\tcellA\tcellB\nMT-ND1\t1\t2\nNDUFS1\t3\t4\n");
    let data = Reader::with_options(
        &p,
        ReaderOptions {
            strict: true,
            force_format: Some(DetectedFormat::DenseTsvCsv),
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.barcodes, vec!["cellA", "cellB"]);
    assert_eq!(data.metadata.gene_symbols, vec!["MT-ND1", "NDUFS1"]);
}

/// strict mode rejects non-finite numeric tokens.
#[test]
fn s5_strict_rejects_non_finite() {
    let p = temp_file("s5", "tsv", "gene\tC1\tC2\nMT\tInf\t1\n");
    let err = Reader::with_options(
        &p,
        ReaderOptions {
            strict: true,
            force_format: Some(DetectedFormat::DenseTsvCsv),
        },
    )
    .read_all()
    .unwrap_err();
    assert_eq!(err.code, kira_scio::ErrorCode::ValidationError);
}

/// a `Reader` should never double-parse — easy to verify via `read_all`
/// success.
#[test]
fn p1_read_all_succeeds_on_single_pass() {
    let d = temp_dir("p1");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1\n2 2 2\n",
    );
    write(&d.join("features.tsv"), "ENSG1\tG1\nENSG2\tG2\n");
    write(&d.join("barcodes.tsv"), "C1\nC2\n");

    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 2);
    assert_eq!(data.matrix.row_idx.len(), 2);
}
