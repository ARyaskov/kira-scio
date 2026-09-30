//! Property tests for the canonical CSC invariants: every matrix built from
//! arbitrary triplets validates, preserves totals, and stays canonical
//! through gene filtering; the MTX writer/reader round trip agrees with the
//! in-memory builder.

mod common;

use std::collections::BTreeMap;

use kira_scio::{Marginals, MatrixStats, Reader, ReaderOptions, SoaCscMatrix};
use proptest::prelude::*;

/// Small integer-valued triplets over a small shape; zeros and negatives
/// included so cancellation and sign handling are exercised.
fn triplets() -> impl Strategy<Value = (usize, usize, Vec<(u32, u32, f32)>)> {
    (1usize..6, 1usize..6).prop_flat_map(|(n_cells, n_genes)| {
        let entry = (
            0..n_cells as u32,
            0..n_genes as u32,
            (-3i32..=5).prop_map(|v| v as f32),
        );
        (
            Just(n_cells),
            Just(n_genes),
            prop::collection::vec(entry, 0..24),
        )
    })
}

/// Reference: dense accumulation of the triplets.
fn reference(triplets: &[(u32, u32, f32)]) -> BTreeMap<(u32, u32), f32> {
    let mut m = BTreeMap::new();
    for &(c, r, v) in triplets {
        *m.entry((c, r)).or_insert(0.0) += v;
    }
    m.retain(|_, v| *v != 0.0);
    m
}

fn matrix_from_reader(dir: &std::path::Path) -> SoaCscMatrix {
    Reader::with_options(
        dir,
        ReaderOptions {
            strict: false,
            ..Default::default()
        },
    )
    .read_matrix()
    .unwrap()
}

/// Writes an MTX triplet (genes x cells) from `(col, row, value)` entries.
fn write_mtx(dir: &common::TestDir, n_cells: usize, n_genes: usize, t: &[(u32, u32, f32)]) {
    let mut body = format!(
        "%%MatrixMarket matrix coordinate real general\n{n_genes} {n_cells} {}\n",
        t.len()
    );
    for &(c, r, v) in t {
        body.push_str(&format!("{} {} {}\n", r + 1, c + 1, v));
    }
    common::write(&dir.join("matrix.mtx"), &body);
    let features: String = (0..n_genes).map(|g| format!("G{g}\tG{g}\n")).collect();
    let barcodes: String = (0..n_cells).map(|c| format!("C{c}\n")).collect();
    common::write(&dir.join("features.tsv"), &features);
    common::write(&dir.join("barcodes.tsv"), &barcodes);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn from_triplets_is_canonical_and_preserves_sums((n_cells, n_genes, t) in triplets()) {
        let (m, merged) = SoaCscMatrix::try_from_triplets(n_cells, n_genes, t.clone()).unwrap();
        prop_assert!(m.validate().is_ok());
        let expected = reference(&t);
        // Same set of stored coordinates and values.
        let mut got = BTreeMap::new();
        for (col, w) in m.col_ptr.windows(2).enumerate() {
            for k in w[0] as usize..w[1] as usize {
                got.insert((col as u32, m.row_idx[k]), m.values[k]);
            }
        }
        prop_assert_eq!(&got, &expected);
        // Merged count = entries that shared a coordinate with an earlier one.
        let distinct: std::collections::BTreeSet<_> = t.iter().map(|e| (e.0, e.1)).collect();
        prop_assert_eq!(merged, t.len() - distinct.len());
        // Whole-matrix statistics agree with the reference.
        let stats = MatrixStats::from_matrix(&m);
        prop_assert_eq!(stats.nnz, expected.len());
        let total: f32 = expected.values().sum();
        prop_assert!((stats.total_counts - f64::from(total)).abs() < 1e-3);
        prop_assert_eq!(stats.has_negative, expected.values().any(|v| *v < 0.0));
        prop_assert!(stats.is_integer);
    }

    #[test]
    fn marginals_sum_to_the_same_total((n_cells, n_genes, t) in triplets()) {
        let (m, _) = SoaCscMatrix::try_from_triplets(n_cells, n_genes, t).unwrap();
        let mg = Marginals::from_matrix(&m);
        prop_assert_eq!(mg.cell_total_counts.len(), n_cells);
        prop_assert_eq!(mg.gene_total_counts.len(), n_genes);
        let by_cell: f64 = mg.cell_total_counts.iter().sum();
        let by_gene: f64 = mg.gene_total_counts.iter().sum();
        prop_assert!((by_cell - by_gene).abs() < 1e-6);
        let nnz_by_cell: u32 = mg.cell_n_features.iter().sum();
        let nnz_by_gene: u32 = mg.gene_n_cells.iter().sum();
        prop_assert_eq!(nnz_by_cell as usize, m.values.len());
        prop_assert_eq!(nnz_by_gene as usize, m.values.len());
        for f in &mg.cell_top_fraction {
            prop_assert!(f.is_finite());
        }
        let all = Marginals::fraction_in_gene_set(&m, &vec![true; n_genes]).unwrap();
        for (frac, total) in all.iter().zip(&mg.cell_total_counts) {
            if *total != 0.0 {
                prop_assert!((frac - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn retain_genes_keeps_canonical_form_and_columns(
        (n_cells, n_genes, t) in triplets(),
        seed in any::<u64>(),
    ) {
        let (m, _) = SoaCscMatrix::try_from_triplets(n_cells, n_genes, t).unwrap();
        let keep: Vec<bool> = (0..n_genes).map(|g| (seed >> (g % 64)) & 1 == 1).collect();
        let r = m.retain_genes(&keep).unwrap();
        prop_assert!(r.validate().is_ok());
        prop_assert_eq!(r.n_cells, n_cells);
        prop_assert_eq!(r.n_genes, keep.iter().filter(|k| **k).count());
        // Kept entries are exactly the original ones whose gene is kept.
        let expected: usize = m.row_idx.iter().filter(|&&g| keep[g as usize]).count();
        prop_assert_eq!(r.values.len(), expected);
        let kept_total: f32 = m
            .row_idx
            .iter()
            .zip(&m.values)
            .filter(|(g, _)| keep[**g as usize])
            .map(|(_, v)| *v)
            .sum();
        let got_total: f32 = r.values.iter().sum();
        prop_assert!((kept_total - got_total).abs() < 1e-3);
    }

    #[test]
    fn mtx_round_trip_matches_in_memory_builder((n_cells, n_genes, t) in triplets()) {
        let dir = common::temp_dir("prop_mtx");
        write_mtx(&dir, n_cells, n_genes, &t);
        let from_file = matrix_from_reader(&dir);
        let (in_memory, _) = SoaCscMatrix::try_from_triplets(n_cells, n_genes, t).unwrap();
        prop_assert_eq!(from_file.col_ptr, in_memory.col_ptr);
        prop_assert_eq!(from_file.row_idx, in_memory.row_idx);
        prop_assert_eq!(from_file.values, in_memory.values);
    }
}
