//! 10x Cell Ranger HDF5 matrices; layouts described in
//! `tests/fixtures/tenx_h5/generate.py`.
#![cfg(feature = "tenx-h5")]

use std::path::PathBuf;

use kira_scio::{
    DetectedFormat, ErrorCode, FeatureTypeFilter, Reader, ReaderOptions, detect_input_format,
};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tenx_h5")
        .join(format!("{name}.h5"))
}

fn assert_reference(data: &kira_scio::CanonicalData) {
    assert_eq!(data.metadata.format, "tenx_h5");
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(
        data.metadata.barcodes,
        vec!["AAACCCAAGAAACACT-1", "AAACCCAAGAAACCAT-1"]
    );
    assert_eq!(
        data.metadata.gene_ids,
        vec!["ENSG1", "ENSG2", "CD3_TotalSeqB"]
    );
    assert_eq!(data.metadata.gene_symbols, vec!["A", "B", "CD3_TotalSeqB"]);
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 2, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 1.0, 3.0]);
    assert!(data.metadata.stats.is_integer);
    assert!(data.metadata.report.is_clean());
    data.matrix.validate().unwrap();
}

#[test]
fn detects_h5_by_extension() {
    assert_eq!(
        detect_input_format(&fixture("v3_filtered")).unwrap(),
        DetectedFormat::TenxH5
    );
}

#[test]
fn cell_ranger_v3_layout_with_feature_types() {
    let data = Reader::new(fixture("v3_filtered")).read_all().unwrap();
    assert_reference(&data);
    assert_eq!(
        data.metadata.feature_types.as_deref(),
        Some(&["Gene Expression", "Gene Expression", "Antibody Capture"].map(String::from)[..])
    );
    assert_eq!(data.metadata.marginals.cell_total_counts, vec![6.0, 3.0]);

    let gex = Reader::with_options(
        fixture("v3_filtered"),
        ReaderOptions {
            strict: true,
            force_format: None,
            feature_types: FeatureTypeFilter::GeneExpression,
            h5ad_source: Default::default(),
            ..Default::default()
        },
    )
    .read_all()
    .unwrap();
    assert_eq!(gex.metadata.gene_ids, vec!["ENSG1", "ENSG2"]);
    assert_eq!(gex.matrix.values, vec![5.0, 3.0]);
    assert_eq!(gex.metadata.report.excluded_features, 1);
}

#[test]
fn cell_ranger_v2_single_genome_layout() {
    let data = Reader::new(fixture("v2_legacy")).read_all().unwrap();
    assert_reference(&data);
    assert_eq!(data.metadata.feature_types, None);
}

#[test]
fn multi_genome_v2_files_are_rejected_clearly() {
    let err = Reader::new(fixture("v2_two_genomes"))
        .read_all()
        .unwrap_err();
    assert_eq!(err.code, ErrorCode::UnsupportedFormat);
    assert!(err.message.contains("GRCh38"), "{}", err.message);
    assert!(err.message.contains("mm10"), "{}", err.message);
}

#[test]
fn read_metadata_and_read_matrix_agree_with_read_all() {
    let md = Reader::new(fixture("v3_filtered")).read_metadata().unwrap();
    let mx = Reader::new(fixture("v3_filtered")).read_matrix().unwrap();
    assert_eq!(md.n_cells, mx.n_cells);
    assert_eq!(md.stats.nnz, mx.values.len());
}

#[test]
fn read_shape_matches_read_all() {
    let r = Reader::new(fixture("v3_filtered"));
    assert_eq!(r.read_shape().unwrap(), (2, 3));
    let r = Reader::with_options(
        fixture("v3_filtered"),
        ReaderOptions {
            strict: true,
            feature_types: FeatureTypeFilter::GeneExpression,
            ..Default::default()
        },
    );
    assert_eq!(r.read_shape().unwrap(), (2, 2));
    assert_eq!(
        Reader::new(fixture("v2_legacy")).read_shape().unwrap(),
        (2, 3)
    );
}

#[test]
fn provenance_records_cell_ranger_dialect_and_kind() {
    use kira_scio::MatrixKind;
    let md = Reader::new(fixture("v3_filtered")).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "cellranger-h5-v3");
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Filtered);
    let md = Reader::new(fixture("v2_legacy")).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "cellranger-h5-v2");
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Unknown);
}
