//! `Reader::read_shape` must agree with `read_all` for every format,
//! including transposed MTX and feature filtering, and `read_metadata`
//! is served from the cache after a full parse.

use std::fs;

use kira_scio::{FeatureTypeFilter, Reader, ReaderOptions};

mod common;
use common::{temp_dir, write};

fn assert_shape_matches(reader: &Reader) {
    let shape = reader.read_shape().unwrap();
    let data = reader.read_all().unwrap();
    assert_eq!(shape, (data.metadata.n_cells, data.metadata.n_genes));
    assert_eq!(shape, (data.matrix.n_cells, data.matrix.n_genes));
}

#[test]
fn mtx_shape_from_header() {
    let d = temp_dir("mtx");
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 1\n3 2 2\n",
    );
    write(&d.join("features.tsv"), "ENSG1\tA\nENSG2\tB\nENSG3\tC\n");
    write(&d.join("barcodes.tsv"), "C1\nC2\n");
    let r = Reader::new(&d);
    assert_eq!(r.read_shape().unwrap(), (2, 3));
    assert_shape_matches(&r);
}

#[test]
fn mtx_shape_honors_transposition_and_feature_filter() {
    let d = temp_dir("mtx_tr");
    // Stored 2 (cells) x 3 (genes); labels identify the orientation.
    write(
        &d.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 3 1\n1 1 4\n",
    );
    write(
        &d.join("features.tsv"),
        "ENSG1\tA\tGene Expression\nENSG2\tB\tGene Expression\nADT\tADT\tAntibody Capture\n",
    );
    write(&d.join("barcodes.tsv"), "C1\nC2\n");
    let r = Reader::new(&d);
    assert_eq!(r.read_shape().unwrap(), (2, 3));
    assert_shape_matches(&r);

    let r = Reader::with_options(
        &d,
        ReaderOptions {
            strict: true,
            feature_types: FeatureTypeFilter::GeneExpression,
            ..Default::default()
        },
    );
    assert_eq!(r.read_shape().unwrap(), (2, 2));
    assert_shape_matches(&r);
}

#[test]
fn dense_shapes_for_both_orientations() {
    let d = temp_dir("dense");
    let gene_major = d.join("g.tsv");
    write(
        &gene_major,
        "#comment\ngene\tC1\tC2\tC3\nG1\t1\t0\t0\nG2\t0\t2\t0\n",
    );
    let r = Reader::new(&gene_major);
    assert_eq!(r.read_shape().unwrap(), (3, 2));
    assert_shape_matches(&r);

    let cell_major = d.join("c.csv");
    write(&cell_major, "Cell_Index,G1,G2\n1,1,0\n2,0,2\n3,0,0\n");
    let r = Reader::new(&cell_major);
    assert_eq!(r.read_shape().unwrap(), (3, 2));
    assert_shape_matches(&r);

    let unlabeled = d.join("u.tsv");
    write(&unlabeled, "\tC1\tC2\nG1\t1\t0\n");
    let r = Reader::new(&unlabeled);
    assert_eq!(r.read_shape().unwrap(), (2, 1));
    assert_shape_matches(&r);
}

#[test]
fn read_metadata_is_served_from_cache_after_a_full_parse() {
    let d = temp_dir("cache");
    let p = d.join("m.tsv");
    write(&p, "gene\tC1\nG1\t1\n");
    let r = Reader::new(&p);
    let first = r.read_metadata().unwrap();
    fs::remove_file(&p).unwrap();
    let second = r.read_metadata().unwrap();
    assert_eq!(first.gene_ids, second.gene_ids);
    assert_eq!(first.stats.nnz, second.stats.nnz);
    // The matrix is never cached, so it needs the file again.
    assert!(r.read_matrix().is_err());
    // A clone starts with an empty cache.
    assert!(r.clone().read_metadata().is_err());
}
