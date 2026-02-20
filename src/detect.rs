use std::path::Path;

use crate::error::{ErrorCode, ScioError, ScioResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedFormat {
    Mtx10x,
    BdRhapsodyWta,
    DenseTsvCsv,
    H5ad,
    Loom,
}

pub fn detect_input_format(path: &Path) -> ScioResult<DetectedFormat> {
    if path.is_dir() {
        if crate::formats::mtx10x::contains_mtx_dataset(path)? {
            return Ok(DetectedFormat::Mtx10x);
        }
        if crate::formats::bd_rhapsody::resolve_bd_input_path(path).is_ok() {
            return Ok(DetectedFormat::BdRhapsodyWta);
        }
        return Err(ScioError::new(
            ErrorCode::UnsupportedFormat,
            format!(
                "failed to auto-detect format in directory {}",
                path.display()
            ),
        ));
    }

    if !path.exists() {
        return Err(ScioError::new(
            ErrorCode::InvalidInputPath,
            format!("input path does not exist: {}", path.display()),
        ));
    }

    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    if name.ends_with(".h5ad") {
        return Ok(DetectedFormat::H5ad);
    }
    if name.ends_with(".loom") {
        return Ok(DetectedFormat::Loom);
    }

    if name == "raw_counts.tsv"
        || name == "raw_counts.tsv.gz"
        || name.ends_with("_raw_counts.tsv")
        || name.ends_with("_raw_counts.tsv.gz")
    {
        return Ok(DetectedFormat::BdRhapsodyWta);
    }

    if name.ends_with(".tsv")
        || name.ends_with(".tsv.gz")
        || name.ends_with(".csv")
        || name.ends_with(".csv.gz")
    {
        return Ok(DetectedFormat::DenseTsvCsv);
    }

    if name.ends_with(".mtx") || name.ends_with(".mtx.gz") {
        return Ok(DetectedFormat::Mtx10x);
    }

    Err(ScioError::new(
        ErrorCode::UnsupportedFormat,
        format!("failed to auto-detect format for {}", path.display()),
    ))
}
