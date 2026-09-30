//! 10x v3 `features.tsv` carries a third `feature_type` column. Multi-modal
//! runs (Feature Barcoding) mix `Gene Expression` with `Antibody Capture`
//! rows whose counts are on a different scale; `FeatureTypeFilter` lets the
//! caller keep one modality.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{FeatureTypeFilter, Reader, ReaderOptions};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_ft_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

/// 3 features x 2 cells: two genes and one antibody with a huge count.
fn multimodal(label: &str) -> PathBuf {
    let d = temp_dir(label);
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 5\n1 1 5\n2 1 7\n3 1 90000\n2 2 1\n3 2 80000\n",
    );
    write(
        &d.join("features.tsv"),
        "ENSG1\tCD3E\tGene Expression\nENSG2\tCD4\tGene Expression\nCD3_TotalSeqB\tCD3_TotalSeqB\tAntibody Capture\n",
    );
    write(&d.join("barcodes.tsv"), "C1-1\nC2-1\n");
    d
}

fn with_filter(filter: FeatureTypeFilter) -> ReaderOptions {
    ReaderOptions {
        strict: true,
        force_format: None,
        feature_types: filter,
        h5ad_source: Default::default(),
        ..Default::default()
    }
}

#[test]
fn feature_types_are_exposed_and_kept_by_default() {
    let data = Reader::new(multimodal("default")).read_all().unwrap();
    assert_eq!(
        data.metadata.feature_types.as_deref(),
        Some(&["Gene Expression", "Gene Expression", "Antibody Capture"].map(String::from)[..])
    );
    assert_eq!(data.metadata.n_genes, 3);
    assert_eq!(data.metadata.stats.max_count, 90000.0);
    assert_eq!(data.metadata.report.excluded_features, 0);
}

#[test]
fn gene_expression_filter_drops_antibody_features() {
    let data = Reader::with_options(
        multimodal("gex"),
        with_filter(FeatureTypeFilter::GeneExpression),
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.n_genes, 2);
    assert_eq!(data.metadata.gene_ids, vec!["ENSG1", "ENSG2"]);
    assert_eq!(data.metadata.gene_symbols, vec!["CD3E", "CD4"]);
    assert_eq!(
        data.metadata.feature_types.as_deref(),
        Some(&["Gene Expression", "Gene Expression"].map(String::from)[..])
    );
    // Matrix renumbered and statistics recomputed on the kept features.
    assert_eq!(data.matrix.n_genes, 2);
    assert_eq!(data.matrix.col_ptr, vec![0, 2, 3]);
    assert_eq!(data.matrix.row_idx, vec![0, 1, 1]);
    assert_eq!(data.matrix.values, vec![5.0, 7.0, 1.0]);
    assert_eq!(data.metadata.stats.nnz, 3);
    assert_eq!(data.metadata.stats.max_count, 7.0);
    assert_eq!(data.metadata.stats.total_counts, 13.0);
    assert_eq!(data.metadata.stats.sparsity, 0.25);
    assert_eq!(data.metadata.report.excluded_features, 1);
    assert!(data.metadata.report.is_clean());
    data.matrix.validate().unwrap();

    // read_metadata and read_matrix see the same filtered view.
    let md = Reader::with_options(
        multimodal("gex_md"),
        with_filter(FeatureTypeFilter::GeneExpression),
    )
    .read_metadata()
    .unwrap();
    assert_eq!(md.n_genes, 2);
}

#[test]
fn only_filter_selects_named_modalities_case_insensitively() {
    let data = Reader::with_options(
        multimodal("only"),
        with_filter(FeatureTypeFilter::Only(vec![
            "antibody capture".to_string(),
        ])),
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.gene_ids, vec!["CD3_TotalSeqB"]);
    assert_eq!(data.matrix.values, vec![90000.0, 80000.0]);
    assert_eq!(data.metadata.report.excluded_features, 2);
}

#[test]
fn filter_matching_nothing_yields_an_empty_gene_axis() {
    let data = Reader::with_options(
        multimodal("none"),
        with_filter(FeatureTypeFilter::Only(vec![
            "CRISPR Guide Capture".to_string(),
        ])),
    )
    .read_all()
    .unwrap();
    assert_eq!(data.metadata.n_genes, 0);
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.matrix.col_ptr, vec![0, 0, 0]);
    assert_eq!(data.metadata.stats.nnz, 0);
    assert_eq!(data.metadata.stats.sparsity, 1.0);
}

#[test]
fn two_column_features_have_no_types_and_filters_are_no_ops() {
    let d = temp_dir("v2");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 1\n2 1 2\n",
    );
    write(&d.join("genes.tsv"), "ENSG1\tA\nENSG2\tB\n");
    write(&d.join("barcodes.tsv"), "C1\n");
    let data = Reader::with_options(&d, with_filter(FeatureTypeFilter::GeneExpression))
        .read_all()
        .unwrap();
    assert_eq!(data.metadata.feature_types, None);
    assert_eq!(data.metadata.n_genes, 2);
    assert_eq!(data.metadata.report.excluded_features, 0);
}

#[test]
fn dense_inputs_carry_no_feature_types() {
    let d = temp_dir("dense");
    let p = d.join("m.tsv");
    write(&p, "gene\tC1\nG1\t1\n");
    let data = Reader::with_options(&p, with_filter(FeatureTypeFilter::GeneExpression))
        .read_all()
        .unwrap();
    assert_eq!(data.metadata.feature_types, None);
    assert_eq!(data.metadata.n_genes, 1);
}
