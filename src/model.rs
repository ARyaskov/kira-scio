#[derive(Debug, Clone)]
pub struct SoaCscMatrix {
    pub n_cells: usize,
    pub n_genes: usize,
    pub col_ptr: Vec<usize>,
    pub row_idx: Vec<usize>,
    pub values: Vec<f32>,
}

impl SoaCscMatrix {
    pub fn validate(&self) -> crate::error::ScioResult<()> {
        if self.col_ptr.len() != self.n_cells + 1 {
            return Err(crate::error::ScioError::new(
                crate::error::ErrorCode::ValidationError,
                "col_ptr length must be n_cells + 1",
            ));
        }
        if self.row_idx.len() != self.values.len() {
            return Err(crate::error::ScioError::new(
                crate::error::ErrorCode::ValidationError,
                "row_idx and values length mismatch",
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
