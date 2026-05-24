//! Canonical in-memory representation.
//!
//! `col_ptr: Vec<u64>` and `row_idx: Vec<u32>` match the on-disk
//! `kira-shared-sc-cache` layout. Cast to `usize` for `sprs::CsMat`.

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
        if let Some(&last) = self.col_ptr.last() {
            if last as usize != self.row_idx.len() {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    "col_ptr tail does not match nnz",
                ));
            }
        }
        let bound = self.n_genes as u32;
        if self.row_idx.iter().any(|&r| r >= bound) {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "row_idx contains an out-of-range gene index",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Default)]
pub struct MatrixStats {
    pub nnz: usize,
    pub total_counts: f64,
    pub min_count: f32,
    pub max_count: f32,
    pub sparsity: f64,
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
}

#[derive(Debug, Clone)]
pub struct CanonicalData {
    pub metadata: InputMetadata,
    pub matrix: SoaCscMatrix,
}
