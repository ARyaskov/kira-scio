//! Canonical in-memory representation.
//!
//! `col_ptr: Vec<u64>` and `row_idx: Vec<u32>` match the on-disk
//! `kira-shared-sc-cache` layout. Cast to `usize` for `sprs::CsMat`.
//!
//! Canonical form: within every column the row indices are strictly
//! increasing (sorted, no duplicate coordinates) and no explicit zeros are
//! stored. Every reader produces this form via [`SoaCscMatrix::from_triplets`]
//! and [`SoaCscMatrix::validate`] rejects anything else.

use crate::error::{ErrorCode, ScioError, ScioResult};

#[derive(Debug, Clone)]
pub struct SoaCscMatrix {
    pub n_cells: usize,
    pub n_genes: usize,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values: Vec<f32>,
}

impl SoaCscMatrix {
    /// Structural validation of the CSC layout and the canonical-form
    /// invariants (monotonic `col_ptr`, strictly increasing rows per column,
    /// in-range row indices).
    pub fn validate(&self) -> ScioResult<()> {
        if self.col_ptr.len() != self.n_cells + 1 {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "col_ptr length must be n_cells + 1",
            ));
        }
        if self.row_idx.len() != self.values.len() {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "row_idx and values length mismatch",
            ));
        }
        if self.col_ptr.first().copied() != Some(0) {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "col_ptr must start at 0",
            ));
        }
        if let Some(&last) = self.col_ptr.last()
            && last as usize != self.row_idx.len()
        {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "col_ptr tail does not match nnz",
            ));
        }
        if self.col_ptr.windows(2).any(|w| w[0] > w[1]) {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "col_ptr must be non-decreasing",
            ));
        }
        let bound = self.n_genes as u32;
        if self.row_idx.iter().any(|&r| r >= bound) {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "row_idx contains an out-of-range gene index",
            ));
        }
        for (col, w) in self.col_ptr.windows(2).enumerate() {
            let rows = &self.row_idx[w[0] as usize..w[1] as usize];
            if rows.windows(2).any(|r| r[0] >= r[1]) {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!("column {col}: row indices must be strictly increasing"),
                ));
            }
        }
        Ok(())
    }

    /// Rejects shapes whose indices cannot be represented by the `u32`
    /// triplet/row-index layout.
    pub fn check_dims(n_cells: usize, n_genes: usize) -> ScioResult<()> {
        if n_cells > u32::MAX as usize || n_genes > u32::MAX as usize {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!("shape {n_cells}x{n_genes} exceeds the u32 index range"),
            ));
        }
        Ok(())
    }

    /// Returns a copy that keeps only the genes (rows) flagged in `keep`,
    /// renumbering the surviving rows in order. `keep.len()` must equal
    /// `n_genes`. Canonical form is preserved.
    pub fn retain_genes(&self, keep: &[bool]) -> ScioResult<Self> {
        if keep.len() != self.n_genes {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!(
                    "retain_genes mask has {} entries, expected {}",
                    keep.len(),
                    self.n_genes
                ),
            ));
        }
        let mut new_index: Vec<u32> = vec![u32::MAX; self.n_genes];
        let mut next = 0u32;
        for (i, &k) in keep.iter().enumerate() {
            if k {
                new_index[i] = next;
                next += 1;
            }
        }
        let mut col_ptr: Vec<u64> = Vec::with_capacity(self.n_cells + 1);
        let mut row_idx: Vec<u32> = Vec::with_capacity(self.row_idx.len());
        let mut values: Vec<f32> = Vec::with_capacity(self.values.len());
        col_ptr.push(0);
        for w in self.col_ptr.windows(2) {
            for k in w[0] as usize..w[1] as usize {
                let mapped = new_index[self.row_idx[k] as usize];
                if mapped != u32::MAX {
                    row_idx.push(mapped);
                    values.push(self.values[k]);
                }
            }
            col_ptr.push(row_idx.len() as u64);
        }
        Ok(Self {
            n_cells: self.n_cells,
            n_genes: next as usize,
            col_ptr,
            row_idx,
            values,
        })
    }

    /// Builds a canonical CSC matrix from `(col, row, value)` triplets.
    ///
    /// Duplicate coordinates are summed (Matrix Market / SciPy convention) and
    /// entries whose value is exactly zero after merging are dropped. Returns
    /// the matrix and the number of merged duplicate entries. Callers must
    /// ensure `col < n_cells` and `row < n_genes`.
    pub(crate) fn from_triplets(
        n_cells: usize,
        n_genes: usize,
        mut triplets: Vec<(u32, u32, f32)>,
    ) -> (Self, usize) {
        // Stable: equal coordinates keep file order, so summing duplicates
        // is deterministic bit-for-bit.
        triplets.sort_by_key(|t| (t.0, t.1));

        let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
        let mut row_idx: Vec<u32> = Vec::with_capacity(triplets.len());
        let mut values: Vec<f32> = Vec::with_capacity(triplets.len());
        let mut merged = 0usize;

        col_ptr.push(0);
        let mut cur_col: u32 = 0;
        let mut i = 0;
        while i < triplets.len() {
            let (col, row, mut val) = triplets[i];
            let mut j = i + 1;
            while j < triplets.len() && triplets[j].0 == col && triplets[j].1 == row {
                val += triplets[j].2;
                j += 1;
            }
            merged += j - i - 1;
            i = j;
            if val == 0.0 {
                continue;
            }
            while cur_col < col {
                col_ptr.push(row_idx.len() as u64);
                cur_col += 1;
            }
            row_idx.push(row);
            values.push(val);
        }
        while col_ptr.len() <= n_cells {
            col_ptr.push(row_idx.len() as u64);
        }

        (
            Self {
                n_cells,
                n_genes,
                col_ptr,
                row_idx,
                values,
            },
            merged,
        )
    }
}

/// Whole-matrix summary computed over stored (non-zero) entries.
#[derive(Debug, Clone, Default)]
pub struct MatrixStats {
    pub nnz: usize,
    pub total_counts: f64,
    /// Minimum over stored entries; 0 when the matrix is empty.
    pub min_count: f32,
    /// Maximum over stored entries; 0 when the matrix is empty.
    pub max_count: f32,
    /// `1 - nnz / (n_cells * n_genes)`; 1 for an empty shape.
    pub sparsity: f64,
    /// Every stored value is a whole number. Raw UMI/read counts are; a
    /// `false` here means the matrix was normalized, scaled or otherwise
    /// transformed upstream and must not be treated as counts.
    pub is_integer: bool,
    /// At least one stored value is negative (scaled/centered data).
    pub has_negative: bool,
}

impl MatrixStats {
    pub fn from_matrix(matrix: &SoaCscMatrix) -> Self {
        let nnz = matrix.values.len();
        let mut total = 0f64;
        let mut min = f32::MAX;
        let mut max = f32::MIN;
        let mut is_integer = true;
        for &v in &matrix.values {
            total += f64::from(v);
            min = min.min(v);
            max = max.max(v);
            is_integer &= v.fract() == 0.0;
        }
        let denom = (matrix.n_cells as u64).saturating_mul(matrix.n_genes as u64);
        let sparsity = if denom == 0 {
            1.0
        } else {
            1.0 - (nnz as f64 / denom as f64)
        };
        Self {
            nnz,
            total_counts: total,
            min_count: if nnz > 0 { min } else { 0.0 },
            max_count: if nnz > 0 { max } else { 0.0 },
            sparsity,
            is_integer,
            has_negative: nnz > 0 && min < 0.0,
        }
    }
}

/// Number of highest-count features used for [`Marginals::cell_top_fraction`].
pub const TOP_FEATURES: usize = 50;

/// Per-cell and per-gene marginals: the first-pass QC covariates every
/// single-cell workflow derives before any biology-specific step (library
/// size, detected features, concentration in top features; per-gene
/// detection). Computed over stored entries of the canonical matrix.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Marginals {
    /// Sum of stored values per cell (library size / count depth).
    pub cell_total_counts: Vec<f64>,
    /// Number of stored (non-zero) features per cell.
    pub cell_n_features: Vec<u32>,
    /// Fraction of each cell's total held by its [`TOP_FEATURES`] largest
    /// entries; 0 for an empty cell. Only meaningful for non-negative data.
    pub cell_top_fraction: Vec<f32>,
    /// Sum of stored values per gene.
    pub gene_total_counts: Vec<f64>,
    /// Number of cells in which each gene has a stored entry.
    pub gene_n_cells: Vec<u32>,
}

impl Marginals {
    pub fn from_matrix(matrix: &SoaCscMatrix) -> Self {
        let n_cells = matrix.n_cells;
        let mut cell_total_counts = Vec::with_capacity(n_cells);
        let mut cell_n_features = Vec::with_capacity(n_cells);
        let mut cell_top_fraction = Vec::with_capacity(n_cells);
        let mut gene_total_counts = vec![0f64; matrix.n_genes];
        let mut gene_n_cells = vec![0u32; matrix.n_genes];
        let mut scratch: Vec<f32> = Vec::new();

        for w in matrix.col_ptr.windows(2) {
            let (start, end) = (w[0] as usize, w[1] as usize);
            let rows = &matrix.row_idx[start..end];
            let vals = &matrix.values[start..end];
            let mut total = 0f64;
            for (&r, &v) in rows.iter().zip(vals) {
                total += f64::from(v);
                gene_total_counts[r as usize] += f64::from(v);
                gene_n_cells[r as usize] += 1;
            }
            cell_total_counts.push(total);
            cell_n_features.push(vals.len() as u32);
            cell_top_fraction.push(top_fraction(vals, total, &mut scratch));
        }

        Self {
            cell_total_counts,
            cell_n_features,
            cell_top_fraction,
            gene_total_counts,
            gene_n_cells,
        }
    }

    /// Fraction of each cell's total held by the genes flagged in `mask`
    /// (for example mitochondrial or ribosomal genes selected by the caller
    /// from the symbols). 0 for an empty cell. `mask.len()` must equal
    /// `n_genes`.
    pub fn fraction_in_gene_set(matrix: &SoaCscMatrix, mask: &[bool]) -> ScioResult<Vec<f32>> {
        if mask.len() != matrix.n_genes {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!(
                    "gene mask has {} entries, expected {}",
                    mask.len(),
                    matrix.n_genes
                ),
            ));
        }
        let mut out = Vec::with_capacity(matrix.n_cells);
        for w in matrix.col_ptr.windows(2) {
            let (start, end) = (w[0] as usize, w[1] as usize);
            let mut total = 0f64;
            let mut in_set = 0f64;
            for (&r, &v) in matrix.row_idx[start..end]
                .iter()
                .zip(&matrix.values[start..end])
            {
                total += f64::from(v);
                if mask[r as usize] {
                    in_set += f64::from(v);
                }
            }
            out.push(if total == 0.0 {
                0.0
            } else {
                (in_set / total) as f32
            });
        }
        Ok(out)
    }
}

fn top_fraction(vals: &[f32], total: f64, scratch: &mut Vec<f32>) -> f32 {
    if total == 0.0 || vals.is_empty() {
        return 0.0;
    }
    if vals.len() <= TOP_FEATURES {
        return 1.0;
    }
    scratch.clear();
    scratch.extend_from_slice(vals);
    // Partition so the TOP_FEATURES largest values sit at the tail.
    let pivot = scratch.len() - TOP_FEATURES;
    scratch.select_nth_unstable_by(pivot, |a, b| a.total_cmp(b));
    let top: f64 = scratch[pivot..].iter().map(|&v| f64::from(v)).sum();
    (top / total) as f32
}

/// A declared-versus-observed count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CountMismatch {
    pub expected: usize,
    pub found: usize,
}

/// Everything a reader had to drop, merge or synthesize to produce the
/// canonical matrix.
///
/// The rules are the same for every format:
///
/// * **Lossy repairs** (`dropped_*`, `entry_count_mismatch`, `relabeled_*`)
///   are errors in strict mode and are only recorded in lenient mode.
/// * **Lossless normalizations** (`merged_duplicates`, `explicit_zeros`,
///   `duplicate_*` labels, `bom_stripped`) are recorded in both modes.
/// * Structural corruption (inconsistent index arrays, unreadable headers)
///   is always an error.
///
/// In strict mode every `dropped_*` counter and every mismatch is therefore
/// zero / `None`; the report is still worth persisting alongside the data
/// for provenance.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct IngestReport {
    /// Entries whose coordinate lies outside the declared shape.
    pub dropped_out_of_range: usize,
    /// Entries whose value is NaN or infinite.
    pub dropped_non_finite: usize,
    /// Explicitly stored zeros in a sparse input; never stored in the
    /// canonical matrix. Always 0 for dense inputs.
    pub explicit_zeros: usize,
    /// Duplicate coordinates that were summed into one entry.
    pub merged_duplicates: usize,
    /// Gene ids seen more than once; each occurrence keeps its own row.
    pub duplicate_gene_ids: Vec<String>,
    /// Barcodes seen more than once; each occurrence keeps its own column.
    pub duplicate_barcodes: Vec<String>,
    /// Gene label vector length disagreed with the matrix and was resized
    /// with synthesized labels.
    pub relabeled_genes: Option<CountMismatch>,
    /// Barcode vector length disagreed with the matrix and was resized with
    /// synthesized labels.
    pub relabeled_barcodes: Option<CountMismatch>,
    /// Matrix Market header entry count disagreed with the entries found.
    pub entry_count_mismatch: Option<CountMismatch>,
    /// A UTF-8 byte order mark was removed from at least one text input.
    pub bom_stripped: bool,
    /// The matrix was stored cells x genes and was transposed to the
    /// canonical genes x cells orientation because the feature and barcode
    /// label counts matched the swapped dimensions.
    pub transposed: bool,
    /// Features removed by the caller's [`crate::FeatureTypeFilter`]. This is
    /// a requested selection, not a repair, so it affects neither
    /// [`Self::is_lossless`] nor [`Self::is_clean`].
    pub excluded_features: usize,
}

impl IngestReport {
    /// True when nothing was dropped or synthesized; label duplicates,
    /// merged coordinates, explicit zeros and a BOM do not count as loss.
    pub fn is_lossless(&self) -> bool {
        self.dropped_out_of_range == 0
            && self.dropped_non_finite == 0
            && self.relabeled_genes.is_none()
            && self.relabeled_barcodes.is_none()
            && self.entry_count_mismatch.is_none()
    }

    /// True when the report has nothing to say at all.
    pub fn is_clean(&self) -> bool {
        self.is_lossless()
            && self.explicit_zeros == 0
            && self.merged_duplicates == 0
            && self.duplicate_gene_ids.is_empty()
            && self.duplicate_barcodes.is_empty()
            && !self.bom_stripped
            && !self.transposed
    }
}

#[derive(Debug, Clone)]
pub struct InputMetadata {
    pub format: String,
    pub n_cells: usize,
    pub n_genes: usize,
    pub gene_ids: Vec<String>,
    pub gene_symbols: Vec<String>,
    pub barcodes: Vec<String>,
    pub stats: MatrixStats,
    /// Per-cell and per-gene marginals of the canonical matrix.
    pub marginals: Marginals,
    /// Per-feature modality as declared by the source (10x `features.tsv`
    /// third column, anndata `var/feature_types`), e.g. `Gene Expression`,
    /// `Antibody Capture`. `None` when the source carries no such column.
    pub feature_types: Option<Vec<String>>,
    /// What the reader tolerated or repaired; see [`IngestReport`].
    pub report: IngestReport,
}

/// Cheap probe of an input's dimensions, read from headers and label files
/// without parsing matrix entries. Feature types are carried so the
/// modality filter can be applied to the gene count.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct ShapeProbe {
    pub n_cells: usize,
    pub n_genes: usize,
    pub feature_types: Option<Vec<String>>,
}

#[derive(Debug, Clone)]
pub struct CanonicalData {
    pub metadata: InputMetadata,
    pub matrix: SoaCscMatrix,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_triplets_sums_duplicates_and_sorts_rows() {
        let (m, merged) = SoaCscMatrix::from_triplets(
            2,
            3,
            vec![(1, 2, 1.0), (0, 1, 2.0), (0, 1, 3.0), (0, 0, 4.0)],
        );
        assert_eq!(merged, 1);
        assert_eq!(m.col_ptr, vec![0, 2, 3]);
        assert_eq!(m.row_idx, vec![0, 1, 2]);
        assert_eq!(m.values, vec![4.0, 5.0, 1.0]);
        m.validate().unwrap();
    }

    #[test]
    fn from_triplets_drops_entries_that_cancel_to_zero() {
        let (m, merged) = SoaCscMatrix::from_triplets(1, 1, vec![(0, 0, 2.0), (0, 0, -2.0)]);
        assert_eq!(merged, 1);
        assert!(m.values.is_empty());
        assert_eq!(m.col_ptr, vec![0, 0]);
        m.validate().unwrap();
    }

    #[test]
    fn from_triplets_handles_empty_and_trailing_empty_columns() {
        let (m, _) = SoaCscMatrix::from_triplets(3, 2, vec![(0, 1, 1.0)]);
        assert_eq!(m.col_ptr, vec![0, 1, 1, 1]);
        m.validate().unwrap();
        let (e, _) = SoaCscMatrix::from_triplets(0, 0, vec![]);
        assert_eq!(e.col_ptr, vec![0]);
        e.validate().unwrap();
    }

    #[test]
    fn validate_rejects_non_canonical_layouts() {
        let base = SoaCscMatrix {
            n_cells: 1,
            n_genes: 3,
            col_ptr: vec![0, 2],
            row_idx: vec![0, 1],
            values: vec![1.0, 1.0],
        };
        base.validate().unwrap();

        let dup = SoaCscMatrix {
            row_idx: vec![1, 1],
            ..base.clone()
        };
        assert_eq!(dup.validate().unwrap_err().code, ErrorCode::ValidationError);

        let unsorted = SoaCscMatrix {
            row_idx: vec![2, 0],
            ..base.clone()
        };
        assert!(unsorted.validate().is_err());

        let bad_ptr = SoaCscMatrix {
            n_cells: 2,
            col_ptr: vec![0, 2, 1],
            row_idx: vec![0, 1],
            values: vec![1.0, 1.0],
            ..base.clone()
        };
        assert!(bad_ptr.validate().is_err());

        let bad_start = SoaCscMatrix {
            col_ptr: vec![1, 2],
            ..base.clone()
        };
        assert!(bad_start.validate().is_err());
    }

    #[test]
    fn retain_genes_renumbers_rows_and_keeps_canonical_form() {
        // 2 cells x 4 genes; keep genes 1 and 3.
        let (m, _) = SoaCscMatrix::from_triplets(
            2,
            4,
            vec![
                (0, 0, 1.0),
                (0, 1, 2.0),
                (0, 3, 3.0),
                (1, 2, 4.0),
                (1, 3, 5.0),
            ],
        );
        let r = m.retain_genes(&[false, true, false, true]).unwrap();
        assert_eq!(r.n_genes, 2);
        assert_eq!(r.n_cells, 2);
        assert_eq!(r.col_ptr, vec![0, 2, 3]);
        assert_eq!(r.row_idx, vec![0, 1, 1]);
        assert_eq!(r.values, vec![2.0, 3.0, 5.0]);
        r.validate().unwrap();

        let none = m.retain_genes(&[false; 4]).unwrap();
        assert_eq!(none.n_genes, 0);
        assert_eq!(none.col_ptr, vec![0, 0, 0]);
        none.validate().unwrap();

        assert!(m.retain_genes(&[true; 3]).is_err());
    }

    #[test]
    fn report_classification() {
        let clean = IngestReport::default();
        assert!(clean.is_clean() && clean.is_lossless());
        let merged = IngestReport {
            merged_duplicates: 2,
            ..Default::default()
        };
        assert!(merged.is_lossless() && !merged.is_clean());
        let lossy = IngestReport {
            dropped_non_finite: 1,
            ..Default::default()
        };
        assert!(!lossy.is_lossless());
    }

    #[test]
    fn check_dims_rejects_u32_overflow() {
        SoaCscMatrix::check_dims(10, 10).unwrap();
        assert!(SoaCscMatrix::check_dims(u32::MAX as usize + 1, 1).is_err());
    }

    #[test]
    fn stats_from_matrix() {
        let (m, _) = SoaCscMatrix::from_triplets(2, 2, vec![(0, 0, 3.0), (1, 1, 1.0)]);
        let s = MatrixStats::from_matrix(&m);
        assert_eq!(s.nnz, 2);
        assert_eq!(s.total_counts, 4.0);
        assert_eq!(s.min_count, 1.0);
        assert_eq!(s.max_count, 3.0);
        assert_eq!(s.sparsity, 0.5);
        assert!(s.is_integer);
        assert!(!s.has_negative);
        let empty = MatrixStats::from_matrix(&SoaCscMatrix::from_triplets(0, 5, vec![]).0);
        assert_eq!(empty.sparsity, 1.0);
        assert_eq!(empty.min_count, 0.0);
        assert!(empty.is_integer);
        assert!(!empty.has_negative);
    }

    #[test]
    fn stats_flag_non_integer_and_negative_values() {
        let (m, _) = SoaCscMatrix::from_triplets(1, 2, vec![(0, 0, 1.5), (0, 1, 2.0)]);
        let s = MatrixStats::from_matrix(&m);
        assert!(!s.is_integer);
        assert!(!s.has_negative);
        let (m, _) = SoaCscMatrix::from_triplets(1, 2, vec![(0, 0, -1.0), (0, 1, 2.0)]);
        let s = MatrixStats::from_matrix(&m);
        assert!(s.is_integer);
        assert!(s.has_negative);
    }

    #[test]
    fn marginals_on_a_small_matrix() {
        // 3 cells x 3 genes; cell 2 is empty.
        let (m, _) = SoaCscMatrix::from_triplets(
            3,
            3,
            vec![(0, 0, 5.0), (0, 2, 1.0), (1, 1, 3.0), (1, 2, 2.0)],
        );
        let mg = Marginals::from_matrix(&m);
        assert_eq!(mg.cell_total_counts, vec![6.0, 5.0, 0.0]);
        assert_eq!(mg.cell_n_features, vec![2, 2, 0]);
        assert_eq!(mg.cell_top_fraction, vec![1.0, 1.0, 0.0]);
        assert_eq!(mg.gene_total_counts, vec![5.0, 3.0, 3.0]);
        assert_eq!(mg.gene_n_cells, vec![1, 1, 2]);

        let mito = Marginals::fraction_in_gene_set(&m, &[false, false, true]).unwrap();
        assert!((mito[0] - 1.0 / 6.0).abs() < 1e-6);
        assert!((mito[1] - 0.4).abs() < 1e-6);
        assert_eq!(mito[2], 0.0);
        assert!(Marginals::fraction_in_gene_set(&m, &[true; 2]).is_err());
    }

    #[test]
    fn top_fraction_uses_the_largest_entries() {
        // One cell with 60 genes: 59 ones and a single 41 -> top 50 = 41 + 49.
        let mut triplets: Vec<(u32, u32, f32)> = (0..59).map(|g| (0, g, 1.0)).collect();
        triplets.push((0, 59, 41.0));
        let (m, _) = SoaCscMatrix::from_triplets(1, 60, triplets);
        let mg = Marginals::from_matrix(&m);
        assert_eq!(mg.cell_total_counts, vec![100.0]);
        assert_eq!(mg.cell_n_features, vec![60]);
        assert!((mg.cell_top_fraction[0] - 0.90).abs() < 1e-6);
    }
}
