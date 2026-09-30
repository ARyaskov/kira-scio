//! Per-cell / per-gene marginals and count-integrity flags exposed on
//! `InputMetadata` for every format, recomputed after feature filtering.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_scio::{FeatureTypeFilter, Marginals, Reader, ReaderOptions};

fn temp_dir(label: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_scio_marg_{label}_{ts}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write(path: &Path, content: &str) {
    let mut f = fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

#[test]
fn mtx_marginals_and_integer_flags() {
    let d = temp_dir("mtx");
    // 3 genes x 3 cells; cell 3 is empty.
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 3 4\n1 1 5\n3 1 1\n2 2 3\n3 2 2\n",
    );
    write(
        &d.join("features.tsv"),
        "ENSG1\tA\nENSG2\tB\nENSG3\tMT-CO1\n",
    );
    write(&d.join("barcodes.tsv"), "C1\nC2\nC3\n");
    let data = Reader::new(&d).read_all().unwrap();
    let m = &data.metadata.marginals;
    assert_eq!(m.cell_total_counts, vec![6.0, 5.0, 0.0]);
    assert_eq!(m.cell_n_features, vec![2, 2, 0]);
    assert_eq!(m.cell_top_fraction, vec![1.0, 1.0, 0.0]);
    assert_eq!(m.gene_total_counts, vec![5.0, 3.0, 3.0]);
    assert_eq!(m.gene_n_cells, vec![1, 1, 2]);
    assert!(data.metadata.stats.is_integer);
    assert!(!data.metadata.stats.has_negative);

    // Caller-defined gene set (here: symbols starting with "MT-").
    let mask: Vec<bool> = data
        .metadata
        .gene_symbols
        .iter()
        .map(|s| s.starts_with("MT-"))
        .collect();
    let mito = Marginals::fraction_in_gene_set(&data.matrix, &mask).unwrap();
    assert!((mito[0] - 1.0 / 6.0).abs() < 1e-6);
    assert!((mito[1] - 0.4).abs() < 1e-6);
    assert_eq!(mito[2], 0.0);
}

#[test]
fn scaled_dense_input_is_flagged_as_non_count_data() {
    let d = temp_dir("scaled");
    let p = d.join("scaled.csv");
    write(&p, "gene,C1,C2\nG1,-1.5,0.25\nG2,2,0\n");
    let data = Reader::new(&p).read_all().unwrap();
    assert!(!data.metadata.stats.is_integer);
    assert!(data.metadata.stats.has_negative);
    assert_eq!(data.metadata.marginals.cell_total_counts, vec![0.5, 0.25]);
    assert_eq!(data.metadata.marginals.gene_n_cells, vec![2, 1]);
}

#[test]
fn marginals_follow_the_feature_filter() {
    let d = temp_dir("filtered");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 5\n2 1 9000\n",
    );
    write(
        &d.join("features.tsv"),
        "ENSG1\tA\tGene Expression\nADT1\tADT1\tAntibody Capture\n",
    );
    write(&d.join("barcodes.tsv"), "C1\n");
    let unfiltered = Reader::new(&d).read_all().unwrap();
    assert_eq!(
        unfiltered.metadata.marginals.cell_total_counts,
        vec![9005.0]
    );
    assert_eq!(unfiltered.metadata.marginals.gene_n_cells, vec![1, 1]);

    let filtered = Reader::with_options(
        &d,
        ReaderOptions {
            strict: true,
            force_format: None,
            feature_types: FeatureTypeFilter::GeneExpression,
        },
    )
    .read_metadata()
    .unwrap();
    assert_eq!(filtered.marginals.cell_total_counts, vec![5.0]);
    assert_eq!(filtered.marginals.cell_n_features, vec![1]);
    assert_eq!(filtered.marginals.gene_n_cells, vec![1]);
    assert_eq!(filtered.marginals.gene_total_counts, vec![5.0]);
}
