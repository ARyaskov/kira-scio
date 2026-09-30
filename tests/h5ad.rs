//! H5AD reader against fixtures written by anndata/h5py; see
//! `tests/fixtures/h5ad/generate.py` for the layouts and the expected matrix.
#![cfg(feature = "h5ad")]

use std::path::PathBuf;

use kira_scio::{DetectedFormat, ErrorCode, Reader, ReaderOptions, detect_input_format};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/h5ad")
        .join(format!("{name}.h5ad"))
}

fn lenient() -> ReaderOptions {
    ReaderOptions {
        strict: false,
        force_format: None,
        feature_types: Default::default(),
        h5ad_source: Default::default(),
        ..Default::default()
    }
}

/// The reference matrix shared by most fixtures (2 cells x 3 genes).
fn assert_reference_matrix(data: &kira_scio::CanonicalData) {
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(data.metadata.barcodes, vec!["AAAC-1", "AAAG-1"]);
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2", "ENSG3"]);
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 2, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 1.0, 3.0]);
    assert_eq!(data.metadata.stats.nnz, 3);
    assert_eq!(data.metadata.stats.total_counts, 9.0);
    data.matrix.validate().unwrap();
}

#[test]
fn detects_h5ad_by_extension() {
    assert_eq!(
        detect_input_format(&fixture("csr")).unwrap(),
        DetectedFormat::H5ad
    );
}

#[test]
fn modern_layouts_all_yield_the_same_canonical_matrix() {
    for name in ["csr", "csc", "dense", "named_index"] {
        let data = Reader::new(fixture(name)).read_all().unwrap();
        assert_reference_matrix(&data);
        assert_eq!(data.metadata.format, "h5ad");
        assert_eq!(data.metadata.gene_symbols, vec!["A", "B", "C"], "{name}");
        assert!(
            data.metadata.report.is_clean(),
            "{name}: {:?}",
            data.metadata.report
        );
    }
}

#[test]
fn legacy_string_datasets_variable_and_fixed_width() {
    for name in ["legacy_vlen", "legacy_fixed"] {
        let data = Reader::new(fixture(name)).read_all().unwrap();
        assert_reference_matrix(&data);
        assert_eq!(data.metadata.gene_symbols, vec!["A", "B", "C"], "{name}");
    }
}

#[test]
fn categorical_symbols_and_int64_data() {
    let data = Reader::new(fixture("categorical_int64"))
        .read_all()
        .unwrap();
    assert_reference_matrix(&data);
    assert_eq!(data.metadata.gene_symbols, vec!["A", "B", "C"]);
}

#[test]
fn messy_csr_is_canonicalized_losslessly() {
    let data = Reader::new(fixture("messy")).read_all().unwrap();
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 2, 1]);
    assert_eq!(data.matrix.values, vec![6.0, 1.0, 3.0]);
    let r = &data.metadata.report;
    assert_eq!(r.explicit_zeros, 1);
    assert_eq!(r.merged_duplicates, 1);
    assert!(r.is_lossless());
    assert_eq!(data.metadata.stats.nnz, 3);
}

#[test]
fn non_finite_values_follow_strict_semantics() {
    let err = Reader::new(fixture("nonfinite")).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ValidationError);

    let data = Reader::with_options(fixture("nonfinite"), lenient())
        .read_all()
        .unwrap();
    assert_eq!(data.metadata.report.dropped_non_finite, 1);
    assert_eq!(data.matrix.col_ptr, vec![0, 1, 2]);
    assert_eq!(data.matrix.row_idx, vec![1, 1]);
    assert_eq!(data.matrix.values, vec![1.0, 2.0]);
}

#[test]
fn out_of_range_indices_follow_strict_semantics() {
    let err = Reader::new(fixture("out_of_range")).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::ValidationError);
    assert!(err.message.contains("out of range"), "{}", err.message);

    let data = Reader::with_options(fixture("out_of_range"), lenient())
        .read_all()
        .unwrap();
    assert_eq!(data.metadata.report.dropped_out_of_range, 1);
    assert_eq!(data.metadata.stats.nnz, 2);
    data.matrix.validate().unwrap();
}

#[test]
fn structural_corruption_is_an_error_in_both_modes() {
    for strict in [true, false] {
        let err = Reader::with_options(
            fixture("short_indptr"),
            ReaderOptions {
                strict,
                force_format: None,
                feature_types: Default::default(),
                h5ad_source: Default::default(),
                ..Default::default()
            },
        )
        .read_all()
        .unwrap_err();
        assert_eq!(err.code, ErrorCode::ParseError, "strict={strict}");
        assert!(err.message.contains("indptr"), "{}", err.message);
    }
}

#[test]
fn x_is_read_even_when_raw_counts_exist() {
    // Documents current behavior: /X (here log-normalized) is what is read;
    // raw/X is not consulted.
    let data = Reader::new(fixture("raw_layer")).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert!(data.matrix.values.iter().all(|v| v.fract() != 0.0));
}

#[test]
fn gzip_compressed_h5ad_is_rejected_with_a_clear_error() {
    let dir = std::env::temp_dir().join(format!(
        "kira_scio_h5ad_gz_{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("x.h5ad.gz");
    std::fs::write(&p, b"not really gzip").unwrap();
    let err = Reader::new(&p).read_all().unwrap_err();
    assert_eq!(err.code, ErrorCode::UnsupportedFormat);
}

#[test]
fn feature_types_column_is_read_and_filterable() {
    use kira_scio::FeatureTypeFilter;
    let data = Reader::new(fixture("feature_types")).read_all().unwrap();
    assert_eq!(
        data.metadata.feature_types.as_deref(),
        Some(&["Gene Expression", "Gene Expression", "Antibody Capture"].map(String::from)[..])
    );
    assert_eq!(data.metadata.stats.max_count, 9000.0);

    let data = Reader::with_options(
        fixture("feature_types"),
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
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2"]);
    assert_eq!(data.matrix.col_ptr, vec![0, 1, 2]);
    assert_eq!(data.matrix.row_idx, vec![0, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 3.0]);
    assert_eq!(data.metadata.stats.max_count, 5.0);
    assert_eq!(data.metadata.report.excluded_features, 1);
}

#[test]
fn matrix_source_selects_x_raw_x_or_a_layer() {
    use kira_scio::H5adSource;
    let opts = |src: H5adSource| ReaderOptions {
        strict: true,
        force_format: None,
        feature_types: Default::default(),
        h5ad_source: src,
        ..Default::default()
    };

    // /X is log-normalized in this fixture.
    let x = Reader::with_options(fixture("raw_layer"), opts(H5adSource::X))
        .read_all()
        .unwrap();
    assert!(!x.metadata.stats.is_integer);

    // raw/X and layers/counts both hold the reference counts.
    let raw = Reader::with_options(fixture("raw_layer"), opts(H5adSource::RawX))
        .read_all()
        .unwrap();
    assert_reference_matrix(&raw);
    assert!(raw.metadata.stats.is_integer);
    assert_eq!(raw.metadata.gene_symbols, vec!["A", "B", "C"]);

    let layer = Reader::with_options(
        fixture("raw_layer"),
        opts(H5adSource::Layer("counts".to_string())),
    )
    .read_all()
    .unwrap();
    assert_reference_matrix(&layer);
    assert!(layer.metadata.stats.is_integer);

    // Absent sources are reported, not silently substituted.
    let err = Reader::with_options(fixture("csr"), opts(H5adSource::RawX))
        .read_all()
        .unwrap_err();
    assert_eq!(err.code, ErrorCode::MissingFile);
    let err = Reader::with_options(
        fixture("raw_layer"),
        opts(H5adSource::Layer("missing".to_string())),
    )
    .read_all()
    .unwrap_err();
    assert_eq!(err.code, ErrorCode::MissingFile);
}

#[test]
fn read_shape_matches_read_all_for_each_source_and_filter() {
    use kira_scio::{FeatureTypeFilter, H5adSource};
    for name in ["csr", "csc", "dense", "legacy_fixed"] {
        let r = Reader::new(fixture(name));
        assert_eq!(r.read_shape().unwrap(), (2, 3), "{name}");
    }
    let r = Reader::with_options(
        fixture("feature_types"),
        ReaderOptions {
            strict: true,
            feature_types: FeatureTypeFilter::GeneExpression,
            ..Default::default()
        },
    );
    assert_eq!(r.read_shape().unwrap(), (2, 2));
    let data = r.read_all().unwrap();
    assert_eq!((data.metadata.n_cells, data.metadata.n_genes), (2, 2));
    let r = Reader::with_options(
        fixture("raw_layer"),
        ReaderOptions {
            strict: true,
            h5ad_source: H5adSource::RawX,
            ..Default::default()
        },
    );
    assert_eq!(r.read_shape().unwrap(), (2, 3));
}

#[test]
fn provenance_records_dialect_and_matrix_source() {
    use kira_scio::{H5adSource, MatrixKind};
    let md = Reader::new(fixture("csr")).read_metadata().unwrap();
    assert_eq!(md.provenance.dialect, "anndata-0.1.0");
    assert_eq!(md.provenance.matrix_source.as_deref(), Some("X"));
    assert_eq!(md.provenance.matrix_kind, MatrixKind::Unknown);
    assert_eq!(md.provenance.source_path, fixture("csr"));
    let md = Reader::with_options(
        fixture("raw_layer"),
        ReaderOptions {
            strict: true,
            h5ad_source: H5adSource::Layer("counts".to_string()),
            ..Default::default()
        },
    )
    .read_metadata()
    .unwrap();
    assert_eq!(
        md.provenance.matrix_source.as_deref(),
        Some("layers/counts")
    );
}
