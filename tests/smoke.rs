use kira_scio::{DetectedFormat, Reader, detect_input_format};

mod common;
use common::{temp_dir, write};

#[test]
fn detect_mtx_prefixed_dot() {
    let d = temp_dir("mtx_dot");
    write(
        &d.join("ABC.matrix.mtx"),
        "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1\n2 2 2\n",
    );
    write(&d.join("ABC.features.tsv"), "ENSG1\tG1\nENSG2\tG2\n");
    write(&d.join("ABC.barcodes.tsv"), "C1\nC2\n");

    assert_eq!(detect_input_format(&d).unwrap(), DetectedFormat::Mtx10x);
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 2);
}

#[test]
fn detect_bd_raw_counts() {
    let d = temp_dir("bd");
    write(
        &d.join("X_raw_counts.tsv"),
        "gene\tC1\tC2\nG1\t1\t0\nG2\t0\t2\n",
    );
    assert_eq!(
        detect_input_format(&d).unwrap(),
        DetectedFormat::BdRhapsodyWta
    );
    let data = Reader::new(&d).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
}

#[test]
fn dense_cell_major_supported() {
    let d = temp_dir("dense_cell_major");
    let p = d.join("in.tsv");
    write(&p, "cell\tG1\tG2\nC1\t1\t0\nC2\t0\t3\n");
    let data = Reader::new(&p).read_all().unwrap();
    assert_eq!(data.metadata.n_cells, 2);
    assert_eq!(data.metadata.n_genes, 2);
}
