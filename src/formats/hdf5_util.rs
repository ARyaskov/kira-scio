//! Helpers shared by the HDF5-based readers (`h5ad`, 10x `.h5`): string
//! datasets in every encoding HDF5 writers use, scalar string attributes,
//! and conversion of compressed-sparse arrays into the canonical CSC matrix
//! under the crate's strict/lenient rules.

use std::path::Path;

use hdf5::types::{FixedAscii, FixedUnicode, TypeDescriptor, VarLenAscii, VarLenUnicode};
use hdf5::{Dataset, Group, Location};

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{IngestReport, SoaCscMatrix};

pub(crate) fn parse_err(msg: impl Into<String>, source: &Path) -> ScioError {
    ScioError::new(ErrorCode::ParseError, msg).with_path(source.to_path_buf())
}

pub(crate) fn hdf5_err(e: hdf5::Error, source: &Path) -> ScioError {
    parse_err(e.to_string(), source)
}

/// Reads a string dataset regardless of whether it is stored as
/// variable-length or fixed-length, ASCII or UTF-8.
pub(crate) fn read_string_dataset(ds: &Dataset, source: &Path) -> ScioResult<Vec<String>> {
    let desc = ds
        .dtype()
        .and_then(|t| t.to_descriptor())
        .map_err(|e| hdf5_err(e, source))?;
    macro_rules! fixed {
        ($ty:ident, $n:expr) => {{
            // HDF5 converts between fixed-length string widths; pick the
            // smallest bucket that fits so buffers stay small.
            match $n {
                0..=32 => ds
                    .read_raw::<$ty<32>>()
                    .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                33..=128 => ds
                    .read_raw::<$ty<128>>()
                    .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                129..=512 => ds
                    .read_raw::<$ty<512>>()
                    .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                513..=2048 => ds
                    .read_raw::<$ty<2048>>()
                    .map(|v| v.iter().map(|s| s.as_str().to_string()).collect()),
                n => {
                    return Err(parse_err(
                        format!(
                            "fixed-length string width {n} exceeds the supported maximum of 2048"
                        ),
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

/// Reads a scalar string attribute (variable- or fixed-length). `None` when
/// the attribute is absent or not a string.
pub(crate) fn read_attr_string(
    loc: &Location,
    name: &str,
    source: &Path,
) -> ScioResult<Option<String>> {
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

/// h5py stores booleans as an int8-backed enum; fall back to raw int8 when
/// the enum conversion is unavailable.
pub(crate) fn read_bool_dataset(ds: &Dataset, source: &Path) -> ScioResult<Vec<bool>> {
    if let Ok(v) = ds.read_raw::<bool>() {
        return Ok(v);
    }
    ds.read_raw::<i8>()
        .map(|v| v.into_iter().map(|b| b != 0).collect())
        .map_err(|e| hdf5_err(e, source))
}

/// Which axis `indptr` runs over in a compressed-sparse group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MajorAxis {
    /// One `indptr` slot per cell; `indices` are gene indices. AnnData
    /// `csr_matrix` (cells x genes) and 10x `.h5` (features x barcodes, CSC).
    Cells,
    /// One `indptr` slot per gene; `indices` are cell indices. AnnData
    /// `csc_matrix`.
    Genes,
}

/// Builds the canonical CSC matrix from a group holding `indptr`, `indices`
/// and `data` datasets (AnnData sparse `X`, 10x `matrix`). `label` names the
/// group in messages, e.g. `/X`.
///
/// Structural inconsistency between the three arrays is always a
/// `ParseError`. Out-of-range indices and non-finite values follow the
/// strict/lenient rules and are counted in `report`; explicit zeros and
/// duplicate coordinates are normalized and counted in both modes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compressed_to_csc(
    group: &Group,
    label: &str,
    major_axis: MajorAxis,
    n_cells: usize,
    n_genes: usize,
    strict: bool,
    source: &Path,
    report: &mut IngestReport,
) -> ScioResult<SoaCscMatrix> {
    SoaCscMatrix::check_dims(n_cells, n_genes).map_err(|e| e.with_path(source.to_path_buf()))?;

    let read_i64 = |name: &str| -> ScioResult<Vec<i64>> {
        group
            .dataset(name)
            .map_err(|_| parse_err(format!("missing {label}/{name}"), source))?
            .read_raw::<i64>()
            .map_err(|e| hdf5_err(e, source))
    };
    let indptr = read_i64("indptr")?;
    let indices = read_i64("indices")?;
    let data: Vec<f32> = group
        .dataset("data")
        .map_err(|_| parse_err(format!("missing {label}/data"), source))?
        .read_raw()
        .map_err(|e| hdf5_err(e, source))?;

    let (major, minor) = match major_axis {
        MajorAxis::Cells => (n_cells, n_genes),
        MajorAxis::Genes => (n_genes, n_cells),
    };

    // Structural integrity: always an error, independent of strictness.
    if indptr.len() != major + 1 {
        return Err(parse_err(
            format!(
                "{label}/indptr has {} entries, expected {}",
                indptr.len(),
                major + 1
            ),
            source,
        ));
    }
    if indptr.first().copied() != Some(0) || indptr.windows(2).any(|w| w[0] > w[1]) {
        return Err(parse_err(
            format!("{label}/indptr must start at 0 and be non-decreasing"),
            source,
        ));
    }
    if indices.len() != data.len() {
        return Err(parse_err(
            format!(
                "{label}/indices ({}) and {label}/data ({}) length mismatch",
                indices.len(),
                data.len()
            ),
            source,
        ));
    }
    if *indptr.last().unwrap_or(&0) as usize != indices.len() {
        return Err(parse_err(
            format!("{label}/indptr tail does not match {label}/indices length"),
            source,
        ));
    }
    if indices.iter().any(|&i| i < 0) {
        return Err(parse_err(
            format!("{label}/indices contains a negative index"),
            source,
        ));
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
                            "{label}/indices[{k}] = {minor_i} is out of range for shape \
                             {n_cells}x{n_genes} (use strict=false to drop)"
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
                        format!("non-finite value in {label}/data[{k}] (use strict=false to drop)"),
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
            let (cell, gene) = match major_axis {
                MajorAxis::Cells => (major_i, minor_i),
                MajorAxis::Genes => (minor_i, major_i),
            };
            triplets.push((cell as u32, gene as u32, val));
        }
    }

    let (matrix, merged) = SoaCscMatrix::from_triplets(n_cells, n_genes, triplets);
    report.merged_duplicates = merged;
    matrix.validate()?;
    Ok(matrix)
}
