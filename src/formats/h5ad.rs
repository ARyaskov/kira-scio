//! H5AD (AnnData) reader.
//!
//! Supported on-disk layouts:
//!
//! * `/X` as a sparse group (`encoding-type` = `csr_matrix` / `csc_matrix`,
//!   or the legacy `h5sparse_format` attribute) or as a dense 2-D dataset.
//!   Index arrays may be int32 or int64; data may be any integer or float
//!   type (converted to `f32` by HDF5).
//! * `obs` / `var` index resolved through the group's `_index` attribute
//!   (anndata >= 0.7), with conventional names as fallback.
//! * String columns stored as variable- or fixed-length string datasets,
//!   as `nullable-string-array` groups (anndata >= 0.10) or as
//!   `categorical` groups (`codes` + `categories`).
//!
//! Strict/lenient semantics follow [`crate::model::IngestReport`]: structural
//! corruption (inconsistent `indptr`/`indices`/`data`) is always an error;
//! out-of-range indices and non-finite values are errors in strict mode and
//! are dropped and counted in lenient mode; explicit zeros and duplicate
//! coordinates are normalized and counted in both modes.

use std::path::Path;

use crate::error::ScioResult;
#[cfg(not(feature = "h5ad"))]
use crate::error::{ErrorCode, ScioError};
#[cfg(not(feature = "h5ad"))]
use crate::model::ShapeProbe;
use crate::model::{InputMetadata, SoaCscMatrix};

use crate::api::H5adSource;

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let (md, _) = read_all(path, strict, &H5adSource::X)?;
    Ok(md)
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    let (_, mx) = read_all(path, strict, &H5adSource::X)?;
    Ok(mx)
}

#[cfg(not(feature = "h5ad"))]
pub(crate) fn read_all(
    path: &Path,
    _strict: bool,
    _source: &H5adSource,
) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "h5ad feature is disabled for this build",
    )
    .with_path(path.to_path_buf()))
}

#[cfg(not(feature = "h5ad"))]
pub(crate) fn read_shape(path: &Path, _source: &H5adSource) -> ScioResult<ShapeProbe> {
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "h5ad feature is disabled for this build",
    )
    .with_path(path.to_path_buf()))
}

#[cfg(feature = "h5ad")]
pub(crate) use imp::{read_all, read_shape};

#[cfg(feature = "h5ad")]
mod imp {
    use std::path::Path;

    use hdf5::{Dataset, File, Group};
    use tracing::warn;

    use crate::api::H5adSource;
    use crate::error::{ErrorCode, ScioError, ScioResult};
    use crate::formats::hdf5_util::{
        MajorAxis, compressed_to_csc, hdf5_err, parse_err, read_attr_string, read_bool_dataset,
        read_string_dataset,
    };
    use crate::model::{IngestReport, InputMetadata, MatrixStats, ShapeProbe, SoaCscMatrix};
    use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol};

    /// Matrix location and the var group labelling its gene axis.
    fn matrix_and_var_paths(which: &H5adSource) -> (String, &'static str) {
        match which {
            H5adSource::X => ("X".to_string(), "var"),
            H5adSource::RawX => ("raw/X".to_string(), "raw/var"),
            H5adSource::Layer(name) => (format!("layers/{name}"), "var"),
        }
    }

    fn open(path: &Path) -> ScioResult<File> {
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
        {
            return Err(ScioError::new(
                ErrorCode::UnsupportedFormat,
                "gzip-compressed .h5ad cannot be opened by HDF5; decompress it first",
            )
            .with_path(path.to_path_buf()));
        }
        File::open(path)
            .map_err(|e| ScioError::new(ErrorCode::Io, e.to_string()).with_path(path.to_path_buf()))
    }

    /// Shape from the matrix metadata plus the feature-type column; no
    /// entries are read.
    pub(crate) fn read_shape(path: &Path, which: &H5adSource) -> ScioResult<ShapeProbe> {
        let file = open(path)?;
        let (matrix_path, var_path) = matrix_and_var_paths(which);
        let (n_cells, n_genes) = if let Ok(group) = file.group(&matrix_path) {
            read_shape_attr(&group, path)?
        } else if let Ok(ds) = file.dataset(&matrix_path) {
            let dims = ds.shape();
            if dims.len() != 2 {
                return Err(parse_err(format!("dense /{matrix_path} must be 2D"), path));
            }
            (dims[0], dims[1])
        } else {
            return Err(ScioError::new(
                ErrorCode::MissingFile,
                format!("requested matrix /{matrix_path} is not present in the file"),
            )
            .with_path(path.to_path_buf()));
        };
        let feature_types =
            read_var_column(&file, var_path, FEATURE_TYPE_COLUMNS, n_genes, false, path)?;
        Ok(ShapeProbe {
            n_cells,
            n_genes,
            feature_types,
        })
    }

    /// Index column names tried when the group carries no `_index` attribute.
    const BARCODE_FALLBACKS: &[&str] = &[
        "_index", "index", "barcode", "barcodes", "cell_id", "cellid",
    ];
    const GENE_ID_FALLBACKS: &[&str] = &["_index", "index", "gene_ids", "gene_id", "feature_id"];
    /// `var` columns that may hold the 10x feature modality.
    const FEATURE_TYPE_COLUMNS: &[&str] = &["feature_types", "feature_type"];
    /// `var` columns that may hold display symbols, in preference order.
    const GENE_SYMBOL_COLUMNS: &[&str] = &[
        "gene_symbols",
        "feature_name",
        "gene_symbol",
        "gene_name",
        "symbol",
    ];

    pub(crate) fn read_all(
        path: &Path,
        strict: bool,
        which: &H5adSource,
    ) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
        let file = open(path)?;
        let mut report = IngestReport::default();
        let (matrix_path, var_path) = matrix_and_var_paths(which);
        if !file.link_exists(&matrix_path) {
            return Err(ScioError::new(
                ErrorCode::MissingFile,
                format!("requested matrix /{matrix_path} is not present in the file"),
            )
            .with_path(path.to_path_buf()));
        }

        let barcodes = read_index(&file, "obs", BARCODE_FALLBACKS, path)?
            .into_iter()
            .enumerate()
            .map(|(i, b)| normalize_barcode(&b, i))
            .collect::<Vec<_>>();
        let gene_raw_ids = read_index(&file, var_path, GENE_ID_FALLBACKS, path)?;
        let gene_symbols_raw = read_var_column(
            &file,
            var_path,
            GENE_SYMBOL_COLUMNS,
            gene_raw_ids.len(),
            strict,
            path,
        )?;
        let feature_types = read_var_column(
            &file,
            var_path,
            FEATURE_TYPE_COLUMNS,
            gene_raw_ids.len(),
            strict,
            path,
        )?;

        let gene_ids = gene_raw_ids
            .iter()
            .enumerate()
            .map(|(i, g)| normalize_gene_id(g, None, i))
            .collect::<Vec<_>>();
        let gene_symbols = match gene_symbols_raw {
            Some(syms) => syms
                .iter()
                .enumerate()
                .map(|(i, s)| normalize_gene_symbol(&gene_raw_ids[i], Some(s), i))
                .collect(),
            None => gene_raw_ids
                .iter()
                .enumerate()
                .map(|(i, g)| normalize_gene_symbol(g, None, i))
                .collect(),
        };

        let matrix = read_matrix_at(&file, &matrix_path, strict, path, &mut report)?;
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
        let stats = MatrixStats::from_matrix(&matrix);
        crate::formats::mtx10x::log_report(path, &report);
        let metadata = InputMetadata {
            format: "h5ad".to_string(),
            n_cells: matrix.n_cells,
            n_genes: matrix.n_genes,
            gene_ids,
            gene_symbols,
            barcodes,
            stats,
            feature_types,
            marginals: Default::default(),
            report,
        };
        Ok((metadata, matrix))
    }

    /// Reads one string column of an anndata dataframe group. Missing values
    /// (nullable mask, negative categorical code) become empty strings so the
    /// caller can synthesize a label.
    fn read_string_column(group: &Group, name: &str, source: &Path) -> ScioResult<Vec<String>> {
        if let Ok(ds) = group.dataset(name) {
            return read_string_dataset(&ds, source);
        }
        let sub = group.group(name).map_err(|_| {
            parse_err(
                format!("column {name} is neither a dataset nor a group"),
                source,
            )
        })?;
        let enc = read_attr_string(&sub, "encoding-type", source)?.unwrap_or_default();
        match enc.as_str() {
            "nullable-string-array" => {
                let values_ds = sub.dataset("values").map_err(|_| {
                    parse_err(format!("column {name}: missing values dataset"), source)
                })?;
                let mut values = read_string_dataset(&values_ds, source)?;
                if let Ok(mask_ds) = sub.dataset("mask") {
                    let mask = read_bool_dataset(&mask_ds, source)?;
                    if mask.len() != values.len() {
                        return Err(parse_err(
                            format!("column {name}: mask/values length mismatch"),
                            source,
                        ));
                    }
                    for (v, masked) in values.iter_mut().zip(mask) {
                        if masked {
                            v.clear();
                        }
                    }
                }
                Ok(values)
            }
            "categorical" => {
                let categories_ds = sub.dataset("categories").map_err(|_| {
                    parse_err(format!("column {name}: missing categories dataset"), source)
                })?;
                let categories = read_string_dataset(&categories_ds, source)?;
                let codes: Vec<i64> = sub
                    .dataset("codes")
                    .map_err(|_| {
                        parse_err(format!("column {name}: missing codes dataset"), source)
                    })?
                    .read_raw()
                    .map_err(|e| hdf5_err(e, source))?;
                codes
                    .into_iter()
                    .map(|c| {
                        if c < 0 {
                            Ok(String::new())
                        } else {
                            categories.get(c as usize).cloned().ok_or_else(|| {
                                parse_err(
                                    format!("column {name}: categorical code {c} out of range"),
                                    source,
                                )
                            })
                        }
                    })
                    .collect()
            }
            other => Err(ScioError::new(
                ErrorCode::UnsupportedFormat,
                format!("column {name}: unsupported encoding-type `{other}`"),
            )
            .with_path(source.to_path_buf())),
        }
    }

    /// Resolves the index column of `obs` or `var`: the `_index` attribute
    /// first, then conventional names.
    fn read_index(
        file: &File,
        group_name: &str,
        fallbacks: &[&str],
        source: &Path,
    ) -> ScioResult<Vec<String>> {
        let group = file.group(group_name).map_err(|_| {
            ScioError::new(
                ErrorCode::UnsupportedFormat,
                format!(
                    "/{group_name} is not a group (legacy compound-dataset layout is not supported)"
                ),
            )
            .with_path(source.to_path_buf())
        })?;
        let mut candidates: Vec<String> = Vec::with_capacity(fallbacks.len() + 1);
        if let Some(name) = read_attr_string(&group, "_index", source)? {
            candidates.push(name);
        }
        for f in fallbacks {
            if !candidates.iter().any(|c| c == f) {
                candidates.push((*f).to_string());
            }
        }
        for cand in &candidates {
            if group.link_exists(cand) {
                return read_string_column(&group, cand, source);
            }
        }
        Err(parse_err(
            format!("missing index column under /{group_name}; tried {candidates:?}"),
            source,
        ))
    }

    /// Reads the first present `var` string column among `candidates`. A
    /// column that exists but cannot be read or has the wrong length is an
    /// error in strict mode and is skipped with a warning otherwise.
    fn read_var_column(
        file: &File,
        var_path: &str,
        candidates: &[&str],
        n_genes: usize,
        strict: bool,
        source: &Path,
    ) -> ScioResult<Option<Vec<String>>> {
        let Ok(var) = file.group(var_path) else {
            return Ok(None);
        };
        for col in candidates {
            if !var.link_exists(col) {
                continue;
            }
            let outcome = read_string_column(&var, col, source).and_then(|values| {
                if values.len() == n_genes {
                    Ok(values)
                } else {
                    Err(parse_err(
                        format!("var/{col} has {} entries, expected {n_genes}", values.len()),
                        source,
                    ))
                }
            });
            match outcome {
                Ok(values) => return Ok(Some(values)),
                Err(err) if strict => return Err(err),
                Err(err) => {
                    warn!(path = %source.display(), column = col, "ignoring var column: {err}");
                }
            }
        }
        Ok(None)
    }

    // ---------------------------------------------------------------- matrix

    fn read_matrix_at(
        file: &File,
        matrix_path: &str,
        strict: bool,
        source: &Path,
        report: &mut IngestReport,
    ) -> ScioResult<SoaCscMatrix> {
        // The matrix may be a Group (sparse CSR/CSC) or a Dataset (dense).
        if let Ok(group) = file.group(matrix_path) {
            return read_sparse_x(&group, strict, source, report);
        }
        if let Ok(dataset) = file.dataset(matrix_path) {
            return read_dense_x(&dataset, strict, source, report);
        }
        Err(parse_err(
            format!("/{matrix_path} is neither a group nor a dataset"),
            source,
        ))
    }

    fn sparse_layout(group: &Group, source: &Path) -> ScioResult<MajorAxis> {
        if let Some(enc) = read_attr_string(group, "encoding-type", source)? {
            return match enc.as_str() {
                "csr_matrix" => Ok(MajorAxis::Cells),
                "csc_matrix" => Ok(MajorAxis::Genes),
                other => Err(ScioError::new(
                    ErrorCode::UnsupportedFormat,
                    format!("unsupported H5AD matrix encoding: {other}"),
                )
                .with_path(source.to_path_buf())),
            };
        }
        // anndata < 0.7 wrote `h5sparse_format` instead of `encoding-type`.
        if let Some(fmt) = read_attr_string(group, "h5sparse_format", source)? {
            return match fmt.as_str() {
                "csr" => Ok(MajorAxis::Cells),
                "csc" => Ok(MajorAxis::Genes),
                other => Err(ScioError::new(
                    ErrorCode::UnsupportedFormat,
                    format!("unsupported h5sparse_format: {other}"),
                )
                .with_path(source.to_path_buf())),
            };
        }
        Ok(MajorAxis::Cells)
    }

    fn read_shape_attr(group: &Group, source: &Path) -> ScioResult<(usize, usize)> {
        let attr = group
            .attr("shape")
            .or_else(|_| group.attr("h5sparse_shape"))
            .map_err(|_| parse_err("missing /X shape attribute", source))?;
        let shape: Vec<i64> = attr.read_raw().map_err(|e| hdf5_err(e, source))?;
        if shape.len() != 2 || shape.iter().any(|&d| d < 0) {
            return Err(parse_err(format!("invalid /X shape {shape:?}"), source));
        }
        Ok((shape[0] as usize, shape[1] as usize))
    }

    fn read_sparse_x(
        group: &Group,
        strict: bool,
        source: &Path,
        report: &mut IngestReport,
    ) -> ScioResult<SoaCscMatrix> {
        let layout = sparse_layout(group, source)?;
        let (n_cells, n_genes) = read_shape_attr(group, source)?;
        compressed_to_csc(
            group, "/X", layout, n_cells, n_genes, strict, source, report,
        )
    }

    fn read_dense_x(
        dataset: &Dataset,
        strict: bool,
        source: &Path,
        report: &mut IngestReport,
    ) -> ScioResult<SoaCscMatrix> {
        let shape = dataset.shape();
        if shape.len() != 2 {
            return Err(parse_err("dense /X must be 2D", source));
        }
        // AnnData dense layout: (cells, genes), stored row-major.
        let n_cells = shape[0];
        let n_genes = shape[1];
        SoaCscMatrix::check_dims(n_cells, n_genes)
            .map_err(|e| e.with_path(source.to_path_buf()))?;
        let array: Vec<f32> = dataset.read_raw().map_err(|e| hdf5_err(e, source))?;
        if array.len() != n_cells * n_genes {
            return Err(parse_err(
                "dense /X element count does not match its shape",
                source,
            ));
        }

        let mut triplets: Vec<(u32, u32, f32)> = Vec::new();
        for cell in 0..n_cells {
            for gene in 0..n_genes {
                let v = array[cell * n_genes + gene];
                if v == 0.0 {
                    continue;
                }
                if !v.is_finite() {
                    if strict {
                        return Err(ScioError::new(
                            ErrorCode::ValidationError,
                            format!("non-finite value in dense /X at ({cell}, {gene}) (use strict=false to drop)"),
                        )
                        .with_path(source.to_path_buf()));
                    }
                    report.dropped_non_finite += 1;
                    continue;
                }
                triplets.push((cell as u32, gene as u32, v));
            }
        }
        let (matrix, _merged) = SoaCscMatrix::from_triplets(n_cells, n_genes, triplets);
        matrix.validate()?;
        Ok(matrix)
    }
}
