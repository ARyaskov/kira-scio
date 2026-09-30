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
use crate::model::{InputMetadata, SoaCscMatrix};

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let (md, _) = read_all(path, strict)?;
    Ok(md)
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    let (_, mx) = read_all(path, strict)?;
    Ok(mx)
}

#[cfg(not(feature = "h5ad"))]
pub(crate) fn read_all(path: &Path, _strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "h5ad feature is disabled for this build",
    )
    .with_path(path.to_path_buf()))
}

#[cfg(feature = "h5ad")]
pub(crate) use imp::read_all;

#[cfg(feature = "h5ad")]
mod imp {
    use std::path::Path;

    use hdf5::types::{FixedAscii, FixedUnicode, TypeDescriptor, VarLenAscii, VarLenUnicode};
    use hdf5::{Dataset, File, Group, Location};
    use tracing::warn;

    use crate::error::{ErrorCode, ScioError, ScioResult};
    use crate::model::{IngestReport, InputMetadata, MatrixStats, SoaCscMatrix};
    use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol};

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

    fn parse_err(msg: impl Into<String>, source: &Path) -> ScioError {
        ScioError::new(ErrorCode::ParseError, msg).with_path(source.to_path_buf())
    }

    fn hdf5_err(e: hdf5::Error, source: &Path) -> ScioError {
        parse_err(e.to_string(), source)
    }

    pub(crate) fn read_all(path: &Path, strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
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
        let file = File::open(path).map_err(|e| {
            ScioError::new(ErrorCode::Io, e.to_string()).with_path(path.to_path_buf())
        })?;
        let mut report = IngestReport::default();

        let barcodes = read_index(&file, "obs", BARCODE_FALLBACKS, path)?
            .into_iter()
            .enumerate()
            .map(|(i, b)| normalize_barcode(&b, i))
            .collect::<Vec<_>>();
        let gene_raw_ids = read_index(&file, "var", GENE_ID_FALLBACKS, path)?;
        let gene_symbols_raw =
            read_var_column(&file, GENE_SYMBOL_COLUMNS, gene_raw_ids.len(), strict, path)?;
        let feature_types = read_var_column(
            &file,
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

        let matrix = read_x_matrix(&file, strict, path, &mut report)?;
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
            report,
        };
        Ok((metadata, matrix))
    }

    // ---------------------------------------------------------------- strings

    /// Reads a string dataset regardless of whether it is stored as
    /// variable-length or fixed-length, ASCII or UTF-8.
    fn read_string_dataset(ds: &Dataset, source: &Path) -> ScioResult<Vec<String>> {
        let desc = ds
            .dtype()
            .and_then(|t| t.to_descriptor())
            .map_err(|e| hdf5_err(e, source))?;
        macro_rules! fixed {
            ($ty:ident, $n:expr) => {{
                // HDF5 converts between fixed-length string widths; pick the
                // smallest bucket that fits so buffers stay small.
                match $n {
                    0..=32 => ds.read_raw::<$ty<32>>().map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                    33..=128 => ds.read_raw::<$ty<128>>().map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                    129..=512 => ds.read_raw::<$ty<512>>().map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                    513..=2048 => ds.read_raw::<$ty<2048>>().map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                    n => {
                        return Err(parse_err(
                            format!("fixed-length string width {n} exceeds the supported maximum of 2048"),
                            source,
                        ));
                    }
                }
            }};
        }
        let out: hdf5::Result<Vec<String>> = match desc {
            TypeDescriptor::VarLenUnicode => ds
                .read_raw::<VarLenUnicode>()
                .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
            TypeDescriptor::VarLenAscii => ds
                .read_raw::<VarLenAscii>()
                .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
            TypeDescriptor::FixedAscii(n) => fixed!(FixedAscii, n),
            TypeDescriptor::FixedUnicode(n) => fixed!(FixedUnicode, n),
            other => {
                return Err(parse_err(
                    format!("expected a string dataset, found {other}"),
                    source,
                ));
            }
        };
        out.map_err(|e| hdf5_err(e, source))
    }

    /// Reads a scalar string attribute (variable- or fixed-length).
    fn read_attr_string(loc: &Location, name: &str, source: &Path) -> ScioResult<Option<String>> {
        let Ok(attr) = loc.attr(name) else {
            return Ok(None);
        };
        let desc = attr
            .dtype()
            .and_then(|t| t.to_descriptor())
            .map_err(|e| hdf5_err(e, source))?;
        let value: hdf5::Result<String> = match desc {
            TypeDescriptor::VarLenUnicode => attr
                .read_scalar::<VarLenUnicode>()
                .map(|s| s.as_str().to_string()),
            TypeDescriptor::VarLenAscii => attr
                .read_scalar::<VarLenAscii>()
                .map(|s| s.as_str().to_string()),
            TypeDescriptor::FixedAscii(n) if n <= 512 => attr
                .read_scalar::<FixedAscii<512>>()
                .map(|s| s.as_str().to_string()),
            TypeDescriptor::FixedUnicode(n) if n <= 512 => attr
                .read_scalar::<FixedUnicode<512>>()
                .map(|s| s.as_str().to_string()),
            _ => return Ok(None),
        };
        value.map(Some).map_err(|e| hdf5_err(e, source))
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

    /// h5py stores booleans as an int8-backed enum; fall back to raw int8
    /// when the enum conversion is unavailable.
    fn read_bool_dataset(ds: &Dataset, source: &Path) -> ScioResult<Vec<bool>> {
        if let Ok(v) = ds.read_raw::<bool>() {
            return Ok(v);
        }
        ds.read_raw::<i8>()
            .map(|v| v.into_iter().map(|b| b != 0).collect())
            .map_err(|e| hdf5_err(e, source))
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
        candidates: &[&str],
        n_genes: usize,
        strict: bool,
        source: &Path,
    ) -> ScioResult<Option<Vec<String>>> {
        let Ok(var) = file.group("var") else {
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

    fn read_x_matrix(
        file: &File,
        strict: bool,
        source: &Path,
        report: &mut IngestReport,
    ) -> ScioResult<SoaCscMatrix> {
        // `/X` may be a Group (sparse CSR/CSC) or a Dataset (dense).
        if let Ok(group) = file.group("X") {
            return read_sparse_x(&group, strict, source, report);
        }
        if let Ok(dataset) = file.dataset("X") {
            return read_dense_x(&dataset, strict, source, report);
        }
        Err(parse_err("missing /X (neither group nor dataset)", source))
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum SparseLayout {
        /// Rows are cells (`indptr` has `n_cells + 1` entries).
        Csr,
        /// Rows are genes (`indptr` has `n_genes + 1` entries).
        Csc,
    }

    fn sparse_layout(group: &Group, source: &Path) -> ScioResult<SparseLayout> {
        if let Some(enc) = read_attr_string(group, "encoding-type", source)? {
            return match enc.as_str() {
                "csr_matrix" => Ok(SparseLayout::Csr),
                "csc_matrix" => Ok(SparseLayout::Csc),
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
                "csr" => Ok(SparseLayout::Csr),
                "csc" => Ok(SparseLayout::Csc),
                other => Err(ScioError::new(
                    ErrorCode::UnsupportedFormat,
                    format!("unsupported h5sparse_format: {other}"),
                )
                .with_path(source.to_path_buf())),
            };
        }
        Ok(SparseLayout::Csr)
    }

    fn read_shape(group: &Group, source: &Path) -> ScioResult<(usize, usize)> {
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
        let (n_cells, n_genes) = read_shape(group, source)?;
        SoaCscMatrix::check_dims(n_cells, n_genes)
            .map_err(|e| e.with_path(source.to_path_buf()))?;

        let read_i64 = |name: &str| -> ScioResult<Vec<i64>> {
            group
                .dataset(name)
                .map_err(|_| parse_err(format!("missing /X/{name}"), source))?
                .read_raw::<i64>()
                .map_err(|e| hdf5_err(e, source))
        };
        let indptr = read_i64("indptr")?;
        let indices = read_i64("indices")?;
        let data: Vec<f32> = group
            .dataset("data")
            .map_err(|_| parse_err("missing /X/data", source))?
            .read_raw()
            .map_err(|e| hdf5_err(e, source))?;

        let (major, minor) = match layout {
            SparseLayout::Csr => (n_cells, n_genes),
            SparseLayout::Csc => (n_genes, n_cells),
        };

        // Structural integrity: always an error, independent of strictness.
        if indptr.len() != major + 1 {
            return Err(parse_err(
                format!(
                    "/X/indptr has {} entries, expected {}",
                    indptr.len(),
                    major + 1
                ),
                source,
            ));
        }
        if indptr.first().copied() != Some(0) || indptr.windows(2).any(|w| w[0] > w[1]) {
            return Err(parse_err(
                "/X/indptr must start at 0 and be non-decreasing",
                source,
            ));
        }
        if indices.len() != data.len() {
            return Err(parse_err(
                format!(
                    "/X/indices ({}) and /X/data ({}) length mismatch",
                    indices.len(),
                    data.len()
                ),
                source,
            ));
        }
        if *indptr.last().unwrap_or(&0) as usize != indices.len() {
            return Err(parse_err(
                "/X/indptr tail does not match /X/indices length",
                source,
            ));
        }
        if indices.iter().any(|&i| i < 0) {
            return Err(parse_err("/X/indices contains a negative index", source));
        }

        let mut triplets: Vec<(u32, u32, f32)> = Vec::with_capacity(data.len());
        for major_i in 0..major {
            let start = indptr[major_i] as usize;
            let end = indptr[major_i + 1] as usize;
            for k in start..end {
                let minor_i = indices[k] as usize;
                let val = data[k];
                if minor_i >= minor {
                    if strict {
                        return Err(ScioError::new(
                            ErrorCode::ValidationError,
                            format!(
                                "/X/indices[{k}] = {minor_i} is out of range for shape {n_cells}x{n_genes} \
                                 (use strict=false to drop)"
                            ),
                        )
                        .with_path(source.to_path_buf()));
                    }
                    report.dropped_out_of_range += 1;
                    continue;
                }
                if !val.is_finite() {
                    if strict {
                        return Err(ScioError::new(
                            ErrorCode::ValidationError,
                            format!("non-finite value in /X/data[{k}] (use strict=false to drop)"),
                        )
                        .with_path(source.to_path_buf()));
                    }
                    report.dropped_non_finite += 1;
                    continue;
                }
                if val == 0.0 {
                    report.explicit_zeros += 1;
                    continue;
                }
                let (cell, gene) = match layout {
                    SparseLayout::Csr => (major_i, minor_i),
                    SparseLayout::Csc => (minor_i, major_i),
                };
                triplets.push((cell as u32, gene as u32, val));
            }
        }

        let (matrix, merged) = SoaCscMatrix::from_triplets(n_cells, n_genes, triplets);
        report.merged_duplicates = merged;
        matrix.validate()?;
        Ok(matrix)
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
