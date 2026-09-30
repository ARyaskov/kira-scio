//! 10x Genomics Cell Ranger HDF5 matrix reader
//! (`filtered_feature_bc_matrix.h5`, `raw_feature_bc_matrix.h5`).
//!
//! Cell Ranger >= 3 stores a single `/matrix` group: `barcodes`, `data`,
//! `indices`, `indptr`, `shape` (`[n_features, n_barcodes]`) and a
//! `features` subgroup with `id`, `name`, `feature_type` and per-tag columns.
//! Cell Ranger 2 stores one group per reference genome with `genes` and
//! `gene_names` instead of `features`. The matrix is CSC over barcodes,
//! i.e. already in the crate's canonical orientation.
//!
//! Requires the `tenx-h5` feature.

use std::path::Path;

use crate::error::ScioResult;
#[cfg(not(feature = "tenx-h5"))]
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

#[cfg(not(feature = "tenx-h5"))]
pub(crate) fn read_all(path: &Path, _strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "tenx-h5 feature is disabled for this build",
    )
    .with_path(path.to_path_buf()))
}

#[cfg(feature = "tenx-h5")]
pub(crate) use imp::read_all;

#[cfg(feature = "tenx-h5")]
mod imp {
    use std::path::Path;

    use hdf5::{File, Group};

    use crate::error::{ErrorCode, ScioError, ScioResult};
    use crate::formats::hdf5_util::{
        MajorAxis, compressed_to_csc, hdf5_err, parse_err, read_string_dataset,
    };
    use crate::formats::mtx10x::{fix_length, log_report};
    use crate::model::{IngestReport, InputMetadata, MatrixStats, SoaCscMatrix};
    use crate::normalize::{
        normalize_barcode, normalize_gene_id, normalize_gene_symbol, synth_barcode,
    };

    pub(crate) fn read_all(path: &Path, strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
        let file = File::open(path).map_err(|e| {
            ScioError::new(ErrorCode::Io, e.to_string()).with_path(path.to_path_buf())
        })?;
        let mut report = IngestReport::default();

        let (group, label) = locate_matrix_group(&file, path)?;
        let strings = |name: &str| -> ScioResult<Vec<String>> {
            let ds = group
                .dataset(name)
                .map_err(|_| parse_err(format!("missing {label}/{name}"), path))?;
            read_string_dataset(&ds, path)
        };

        let shape: Vec<i64> = group
            .dataset("shape")
            .map_err(|_| parse_err(format!("missing {label}/shape"), path))?
            .read_raw()
            .map_err(|e| hdf5_err(e, path))?;
        if shape.len() != 2 || shape.iter().any(|&d| d < 0) {
            return Err(parse_err(format!("invalid {label}/shape {shape:?}"), path));
        }
        // Cell Ranger shape is [n_features, n_barcodes].
        let n_genes = shape[0] as usize;
        let n_cells = shape[1] as usize;

        let mut barcodes: Vec<String> = strings("barcodes")?
            .into_iter()
            .enumerate()
            .map(|(i, b)| normalize_barcode(&b, i))
            .collect();

        let (raw_ids, raw_names, feature_types) = if group.link_exists("features") {
            let features = group
                .group("features")
                .map_err(|_| parse_err(format!("{label}/features is not a group"), path))?;
            let read = |name: &str| -> ScioResult<Vec<String>> {
                let ds = features
                    .dataset(name)
                    .map_err(|_| parse_err(format!("missing {label}/features/{name}"), path))?;
                read_string_dataset(&ds, path)
            };
            let ids = read("id")?;
            let names = if features.link_exists("name") {
                read("name")?
            } else {
                ids.clone()
            };
            let types = if features.link_exists("feature_type") {
                Some(read("feature_type")?)
            } else {
                None
            };
            (ids, names, types)
        } else {
            // Cell Ranger 2 layout.
            let ids = strings("genes")?;
            let names = if group.link_exists("gene_names") {
                strings("gene_names")?
            } else {
                ids.clone()
            };
            (ids, names, None)
        };
        if raw_names.len() != raw_ids.len() {
            return Err(parse_err(
                format!(
                    "{label}: {} feature ids but {} feature names",
                    raw_ids.len(),
                    raw_names.len()
                ),
                path,
            ));
        }
        let mut gene_ids: Vec<String> = raw_ids
            .iter()
            .zip(&raw_names)
            .enumerate()
            .map(|(i, (id, name))| normalize_gene_id(id, Some(name), i))
            .collect();
        let mut gene_symbols: Vec<String> = raw_ids
            .iter()
            .zip(&raw_names)
            .enumerate()
            .map(|(i, (id, name))| normalize_gene_symbol(id, Some(name), i))
            .collect();
        let feature_types = feature_types.filter(|t| t.len() == n_genes);

        report.relabeled_genes = fix_length(&mut gene_ids, n_genes, strict, "gene", path, |i| {
            normalize_gene_id("", None, i)
        })?;
        fix_length(
            &mut gene_symbols,
            n_genes,
            strict,
            "gene_symbol",
            path,
            |i| normalize_gene_symbol("", None, i),
        )?;
        report.relabeled_barcodes = fix_length(
            &mut barcodes,
            n_cells,
            strict,
            "barcode",
            path,
            synth_barcode,
        )?;

        let matrix = compressed_to_csc(
            &group,
            &label,
            MajorAxis::Cells,
            n_cells,
            n_genes,
            strict,
            path,
            &mut report,
        )?;
        let stats = MatrixStats::from_matrix(&matrix);
        log_report(path, &report);

        let metadata = InputMetadata {
            format: "tenx_h5".to_string(),
            n_cells,
            n_genes,
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

    /// `/matrix` for Cell Ranger >= 3; otherwise the single per-genome group
    /// of the Cell Ranger 2 layout.
    fn locate_matrix_group(file: &File, source: &Path) -> ScioResult<(Group, String)> {
        if file.link_exists("matrix") {
            let g = file
                .group("matrix")
                .map_err(|_| parse_err("/matrix is not a group", source))?;
            return Ok((g, "/matrix".to_string()));
        }
        let mut genome_groups: Vec<String> = file
            .member_names()
            .map_err(|e| hdf5_err(e, source))?
            .into_iter()
            .filter(|name| {
                file.group(name)
                    .map(|g| g.link_exists("indptr") && g.link_exists("barcodes"))
                    .unwrap_or(false)
            })
            .collect();
        genome_groups.sort();
        match genome_groups.len() {
            1 => {
                let name = genome_groups.remove(0);
                let g = file.group(&name).map_err(|e| hdf5_err(e, source))?;
                Ok((g, format!("/{name}")))
            }
            0 => Err(parse_err(
                "no /matrix group and no Cell Ranger 2 genome group found",
                source,
            )),
            _ => Err(ScioError::new(
                ErrorCode::UnsupportedFormat,
                format!(
                    "multi-genome Cell Ranger 2 matrix ({}) is not supported; \
                     re-run with a single reference or export MEX",
                    genome_groups.join(", ")
                ),
            )
            .with_path(source.to_path_buf())),
        }
    }
}
