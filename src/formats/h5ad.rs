use std::path::Path;

use crate::error::{ErrorCode, ScioError, ScioResult};
#[cfg(feature = "h5ad")]
use crate::model::MatrixStats;
use crate::model::{InputMetadata, SoaCscMatrix};
#[cfg(feature = "h5ad")]
use crate::normalize::{normalize_barcode, normalize_gene_symbol};

#[cfg(feature = "h5ad")]
use hdf5::File;

pub fn read_metadata(path: &Path, _strict: bool) -> ScioResult<InputMetadata> {
    #[cfg(feature = "h5ad")]
    {
        let file = File::open(path).map_err(|e| ScioError::new(ErrorCode::Io, e.to_string()))?;
        let barcodes = read_strings(&file, "obs/_index")?
            .into_iter()
            .enumerate()
            .map(|(i, b)| normalize_barcode(&b, i))
            .collect::<Vec<_>>();
        let genes_raw = read_strings(&file, "var/_index")?;
        let genes = genes_raw
            .iter()
            .enumerate()
            .map(|(i, g)| normalize_gene_symbol(g, Some(g), i))
            .collect::<Vec<_>>();

        let matrix = read_matrix(path, true)?;
        let stats = matrix_stats(&matrix);

        return Ok(InputMetadata {
            format: "h5ad".to_string(),
            n_cells: barcodes.len(),
            n_genes: genes.len(),
            gene_ids: genes.clone(),
            gene_symbols: genes,
            barcodes,
            stats,
        });
    }

    #[cfg(not(feature = "h5ad"))]
    {
        let _ = path;
        Err(ScioError::new(
            ErrorCode::FeatureDisabled,
            "h5ad feature is disabled for this build",
        ))
    }
}

pub fn read_matrix(path: &Path, _strict: bool) -> ScioResult<SoaCscMatrix> {
    #[cfg(feature = "h5ad")]
    {
        let file = File::open(path).map_err(|e| ScioError::new(ErrorCode::Io, e.to_string()))?;
        let x = file
            .group("X")
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X group in h5ad"))?;

        // AnnData sparse CSR/CSC canonical datasets.
        let shape: Vec<usize> = x
            .attr("shape")
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X shape attr"))?
            .read_raw()
            .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;

        if shape.len() != 2 {
            return Err(ScioError::new(ErrorCode::ParseError, "invalid /X shape"));
        }
        let n_cells = shape[0];
        let n_genes = shape[1];

        let indptr: Vec<usize> = x
            .dataset("indptr")
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/indptr"))?
            .read_raw()
            .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
        let indices: Vec<usize> = x
            .dataset("indices")
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/indices"))?
            .read_raw()
            .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;
        let data: Vec<f32> = x
            .dataset("data")
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "missing /X/data"))?
            .read_raw()
            .map_err(|e| ScioError::new(ErrorCode::ParseError, e.to_string()))?;

        let enc =
            read_attr_string(&x, "encoding-type").unwrap_or_else(|_| "csr_matrix".to_string());
        if enc == "csc_matrix" {
            let matrix = SoaCscMatrix {
                n_cells,
                n_genes,
                col_ptr: indptr,
                row_idx: indices,
                values: data,
            };
            matrix.validate()?;
            return Ok(matrix);
        }

        if enc == "csr_matrix" {
            // CSR (cell-major rows) -> CSC conversion.
            let mut cols = vec![Vec::<(usize, f32)>::new(); n_cells];
            for row in 0..n_cells {
                let start = indptr[row];
                let end = indptr[row + 1];
                for idx in start..end {
                    let gene = indices[idx];
                    let val = data[idx];
                    if gene >= n_genes {
                        continue;
                    }
                    cols[row].push((gene, val));
                }
            }

            let mut col_ptr = Vec::with_capacity(n_cells + 1);
            let mut row_idx = Vec::new();
            let mut values = Vec::new();
            col_ptr.push(0);
            for col in &mut cols {
                col.sort_by_key(|(r, _)| *r);
                for (r, v) in col.iter() {
                    row_idx.push(*r);
                    values.push(*v);
                }
                col_ptr.push(row_idx.len());
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

        return Err(ScioError::new(
            ErrorCode::UnsupportedFormat,
            format!("unsupported H5AD matrix encoding: {enc}"),
        ));
    }

    #[cfg(not(feature = "h5ad"))]
    {
        let _ = path;
        Err(ScioError::new(
            ErrorCode::FeatureDisabled,
            "h5ad feature is disabled for this build",
        ))
    }
}

#[cfg(feature = "h5ad")]
fn read_strings(file: &File, path: &str) -> ScioResult<Vec<String>> {
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
    let sparsity = 1.0 - (nnz as f64 / ((matrix.n_cells * matrix.n_genes) as f64));
    MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min } else { 0.0 },
        max_count: if nnz > 0 { max } else { 0.0 },
        sparsity,
    }
}
