//! H5AD (AnnData) reader.
//!
//! Sparse `X` (CSR or CSC) or dense `X` (Dataset). `obs/var` index name
//! falls back across `_index`, `index`, `barcode`, `cell_id`, etc. Inf/NaN
//! in `X/data` errors under strict, otherwise warns and drops.

use std::path::Path;

#[cfg(feature = "h5ad")]
use tracing::warn;

use crate::error::{ErrorCode, ScioError, ScioResult};
#[cfg(feature = "h5ad")]
use crate::model::MatrixStats;
use crate::model::{InputMetadata, SoaCscMatrix};
#[cfg(feature = "h5ad")]
use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol};

#[cfg(feature = "h5ad")]
use hdf5::File;

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    #[cfg(feature = "h5ad")]
    {
        let (md, _) = read_all(path, strict)?;
        return Ok(md);
    }
    #[cfg(not(feature = "h5ad"))]
    {
        let _ = path;
        let _ = strict;
        Err(ScioError::new(
            ErrorCode::FeatureDisabled,
            "h5ad feature is disabled for this build",
        )
        .with_path(path.to_path_buf()))
    }
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    #[cfg(feature = "h5ad")]
    {
        let (_, mx) = read_all(path, strict)?;
        return Ok(mx);
    }
    #[cfg(not(feature = "h5ad"))]
    {
        let _ = path;
        let _ = strict;
        Err(ScioError::new(
            ErrorCode::FeatureDisabled,
            "h5ad feature is disabled for this build",
        )
        .with_path(path.to_path_buf()))
    }
}

#[cfg(feature = "h5ad")]
pub(crate) fn read_all(
    path: &Path,
    strict: bool,
) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    let file = File::open(path).map_err(|e| {
        ScioError::new(ErrorCode::Io, e.to_string())
            .with_path(path.to_path_buf())
    })?;

    let barcodes = read_index_strings(&file, "obs", BARCODE_FALLBACKS)?
        .into_iter()
        .enumerate()
        .map(|(i, b)| normalize_barcode(&b, i))
        .collect::<Vec<_>>();
    let gene_raw_ids = read_index_strings(&file, "var", GENE_ID_FALLBACKS)?;
    let gene_symbols_raw = read_optional_strings(&file, "var/gene_symbols")
        .or_else(|_| read_optional_strings(&file, "var/feature_name"))
        .ok();

    let gene_ids = gene_raw_ids
        .iter()
        .enumerate()
        .map(|(i, g)| normalize_gene_id(g, None, i))
        .collect::<Vec<_>>();
    let gene_symbols = gene_symbols_raw
        .as_ref()
        .map(|syms| {
            syms.iter()
                .enumerate()
                .map(|(i, s)| normalize_gene_symbol(&gene_raw_ids[i], Some(s), i))
                .collect::<Vec<_>>()
        })
        .unwrap_or_else(|| {
            gene_raw_ids
                .iter()
                .enumerate()
                .map(|(i, g)| normalize_gene_symbol(g, None, i))
                .collect()
        });

    let matrix = read_x_matrix(&file, strict, path)?;
    if matrix.n_cells != barcodes.len() || matrix.n_genes != gene_ids.len() {
        return Err(ScioError::new(
            ErrorCode::DimensionMismatch,
            format!(
                "h5ad metadata/matrix mismatch: barcodes={} genes={} matrix={}x{}",
                barcodes.len(),
                gene_ids.len(),
                matrix.n_cells,
                matrix.n_genes
            ),
        )
        .with_path(path.to_path_buf()));
    }
    let stats = matrix_stats(&matrix);
    let metadata = InputMetadata {
        format: "h5ad".to_string(),
        n_cells: matrix.n_cells,
        n_genes: matrix.n_genes,
        gene_ids,
        gene_symbols,
        barcodes,
        stats,
    };
    Ok((metadata, matrix))
}

#[cfg(not(feature = "h5ad"))]
pub(crate) fn read_all(
    path: &Path,
    _strict: bool,
) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "h5ad feature is disabled for this build",
    )
    .with_path(path.to_path_buf()))
}

#[cfg(feature = "h5ad")]
const BARCODE_FALLBACKS: &[&str] = &[
    "_index", "index", "barcode", "barcodes", "cell_id", "cellid",
];
#[cfg(feature = "h5ad")]
const GENE_ID_FALLBACKS: &[&str] = &[
    "_index", "index", "gene_ids", "gene_id", "feature_id",
];

#[cfg(feature = "h5ad")]
fn read_index_strings(
    file: &File,
    group: &str,
    candidates: &[&str],
) -> ScioResult<Vec<String>> {
    use hdf5::types::VarLenUnicode;

    for cand in candidates {
        let dataset_path = format!("{group}/{cand}");
        if let Ok(ds) = file.dataset(&dataset_path) {
            let data: Vec<VarLenUnicode> = ds
                .read_raw()
                .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
            return Ok(data.into_iter().map(|v| v.to_string()).collect());
        }
    }
    Err(ScioError::new(
        ErrorCode::ParseError,
        format!(
            "missing index dataset under {group}; tried {:?}",
            candidates
        ),
    ))
}

#[cfg(feature = "h5ad")]
fn read_optional_strings(file: &File, path: &str) -> ScioResult<Vec<String>> {
    use hdf5::types::VarLenUnicode;
    let ds = file
        .dataset(path)
        .map_err(|_| ScioError::new(ErrorCode::ParseError, format!("missing dataset: {path}")))?;
    let data: Vec<VarLenUnicode> = ds
        .read_raw()
        .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
    Ok(data.into_iter().map(|v| v.to_string()).collect())
}

#[cfg(feature = "h5ad")]
fn read_x_matrix(file: &File, strict: bool, source: &Path) -> ScioResult<SoaCscMatrix> {
    // `/X` may be a Group (sparse CSR/CSC) or a Dataset (dense f32).
    if let Ok(group) = file.group("X") {
        return read_sparse_x(&group, strict, source);
    }
    if let Ok(dataset) = file.dataset("X") {
        return read_dense_x(&dataset, strict, source);
    }
    Err(ScioError::new(
        ErrorCode::ParseError,
        "missing /X (neither group nor dataset)",
    )
    .with_path(source.to_path_buf()))
}

#[cfg(feature = "h5ad")]
fn read_sparse_x(group: &hdf5::Group, strict: bool, source: &Path) -> ScioResult<SoaCscMatrix> {
    let shape: Vec<u64> = group
        .attr("shape")
        .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X shape attr"))?
        .read_raw()
        .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
    if shape.len() != 2 {
        return Err(ScioError::new(ErrorCode::ParseError, "invalid /X shape"));
    }
    let n_cells = shape[0] as usize;
    let n_genes = shape[1] as usize;

    // indptr/indices may be int32 or int64 depending on the writer.
    let indptr_i64: Vec<i64> = group
        .dataset("indptr")
        .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/indptr"))?
        .read_raw::<i64>()
        .or_else(|_| {
            group
                .dataset("indptr")?
                .read_raw::<i32>()
                .map(|v| v.into_iter().map(|x| x as i64).collect())
        })
        .map_err(|e: hdf5::Error| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
    let indices_i64: Vec<i64> = group
        .dataset("indices")
        .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/indices"))?
        .read_raw::<i64>()
        .or_else(|_| {
            group
                .dataset("indices")?
                .read_raw::<i32>()
                .map(|v| v.into_iter().map(|x| x as i64).collect())
        })
        .map_err(|e: hdf5::Error| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
    let data: Vec<f32> = group
        .dataset("data")
        .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/data"))?
        .read_raw()
        .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;

    let enc = read_attr_string(group, "encoding-type").unwrap_or_else(|_| "csr_matrix".to_string());

    if enc == "csc_matrix" {
        let col_ptr: Vec<u64> = indptr_i64.into_iter().map(|v| v as u64).collect();
        let mut row_idx: Vec<u32> = Vec::with_capacity(indices_i64.len());
        let mut values: Vec<f32> = Vec::with_capacity(data.len());
        check_and_collect(
            indices_i64.iter().copied().zip(data.iter().copied()),
            n_genes,
            strict,
            source,
            &mut row_idx,
            &mut values,
        )?;
        let matrix = SoaCscMatrix {
            n_cells,
            n_genes,
            col_ptr,
            row_idx,
            values,
        };
        matrix.validate()?;
        return Ok(matrix);
    }

    if enc == "csr_matrix" {
        // CSR (per-cell rows) → CSC conversion via triplet sort.
        let mut triplets: Vec<(u32, u32, f32)> = Vec::with_capacity(indices_i64.len());
        let mut saw_nonfinite = false;
        for row in 0..n_cells {
            let start = indptr_i64[row] as usize;
            let end = indptr_i64[row + 1] as usize;
            for idx in start..end {
                let gene = indices_i64[idx] as usize;
                let val = data[idx];
                if gene >= n_genes {
                    continue;
                }
                if !val.is_finite() {
                    saw_nonfinite = true;
                    if strict {
                        return Err(ScioError::new(
                            ErrorCode::ValidationError,
                            "non-finite value in /X/data",
                        )
                        .with_path(source.to_path_buf()));
                    }
                    continue;
                }
                triplets.push((row as u32, gene as u32, val));
            }
        }
        if saw_nonfinite {
            warn!(path = %source.display(), "h5ad /X/data: non-finite dropped");
        }
        triplets.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));

        let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
        let mut row_idx: Vec<u32> = Vec::with_capacity(triplets.len());
        let mut values: Vec<f32> = Vec::with_capacity(triplets.len());
        col_ptr.push(0);
        let mut cur_col: u32 = 0;
        for (col, row, val) in triplets {
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

        let matrix = SoaCscMatrix {
            n_cells,
            n_genes,
            col_ptr,
            row_idx,
            values,
        };
        matrix.validate()?;
        return Ok(matrix);
    }

    Err(ScioError::new(
        ErrorCode::UnsupportedFormat,
        format!("unsupported H5AD matrix encoding: {enc}"),
    )
    .with_path(source.to_path_buf()))
}

#[cfg(feature = "h5ad")]
fn check_and_collect(
    pairs: impl Iterator<Item = (i64, f32)>,
    n_genes: usize,
    strict: bool,
    source: &Path,
    row_idx: &mut Vec<u32>,
    values: &mut Vec<f32>,
) -> ScioResult<()> {
    let mut saw_nonfinite = false;
    for (gene, val) in pairs {
        if (gene as usize) >= n_genes {
            continue;
        }
        if !val.is_finite() {
            saw_nonfinite = true;
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    "non-finite value in /X/data",
                )
                .with_path(source.to_path_buf()));
            }
            continue;
        }
        row_idx.push(gene as u32);
        values.push(val);
    }
    if saw_nonfinite {
        warn!(path = %source.display(), "h5ad /X/data: non-finite dropped");
    }
    Ok(())
}

#[cfg(feature = "h5ad")]
fn read_dense_x(
    dataset: &hdf5::Dataset,
    strict: bool,
    source: &Path,
) -> ScioResult<SoaCscMatrix> {
    use ndarray::Ix2;
    let shape = dataset.shape();
    if shape.len() != 2 {
        return Err(ScioError::new(ErrorCode::ParseError, "dense /X must be 2D")
            .with_path(source.to_path_buf()));
    }
    // AnnData dense layout: (cells, genes).
    let n_cells = shape[0];
    let n_genes = shape[1];
    let array = dataset
        .read::<f32, Ix2>()
        .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;

    let mut triplets: Vec<(u32, u32, f32)> = Vec::new();
    let mut saw_nonfinite = false;
    for cell in 0..n_cells {
        for gene in 0..n_genes {
            let v = array[(cell, gene)];
            if v == 0.0 {
                continue;
            }
            if !v.is_finite() {
                saw_nonfinite = true;
                if strict {
                    return Err(ScioError::new(
                        ErrorCode::ValidationError,
                        "non-finite value in dense /X",
                    )
                    .with_path(source.to_path_buf()));
                }
                continue;
            }
            triplets.push((cell as u32, gene as u32, v));
        }
    }
    if saw_nonfinite {
        warn!(path = %source.display(), "dense h5ad /X: non-finite dropped");
    }
    triplets.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));

    let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
    let mut row_idx: Vec<u32> = Vec::with_capacity(triplets.len());
    let mut values: Vec<f32> = Vec::with_capacity(triplets.len());
    col_ptr.push(0);
    let mut cur_col: u32 = 0;
    for (col, row, val) in triplets {
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

    let matrix = SoaCscMatrix {
        n_cells,
        n_genes,
        col_ptr,
        row_idx,
        values,
    };
    matrix.validate()?;
    Ok(matrix)
}

#[cfg(feature = "h5ad")]
fn read_attr_string(group: &hdf5::Group, name: &str) -> ScioResult<String> {
    use hdf5::types::VarLenUnicode;
    let attr = group
        .attr(name)
        .map_err(|_| ScioError::new(ErrorCode::ParseError, format!("missing attr: {name}")))?;
    let value: VarLenUnicode = attr
        .read_scalar()
        .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
    Ok(value.to_string())
}

#[cfg(feature = "h5ad")]
fn matrix_stats(matrix: &SoaCscMatrix) -> MatrixStats {
    let nnz = matrix.values.len();
    let mut total = 0f64;
    let mut min = f32::MAX;
    let mut max = f32::MIN;
    for v in &matrix.values {
        total += *v as f64;
        min = min.min(*v);
        max = max.max(*v);
    }
    let denom = matrix.n_cells.saturating_mul(matrix.n_genes);
    let sparsity = if denom == 0 {
        1.0
    } else {
        1.0 - (nnz as f64 / denom as f64)
    };
    MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min } else { 0.0 },
        max_count: if nnz > 0 { max } else { 0.0 },
        sparsity,
    }
}
