use std::path::{Path, PathBuf};

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, SoaCscMatrix};

pub fn resolve_bd_input_path(path: &Path) -> ScioResult<PathBuf> {
    if path.is_file() {
        return Ok(path.to_path_buf());
    }
    if !path.is_dir() {
        return Err(ScioError::new(
            ErrorCode::InvalidInputPath,
            format!("invalid input path: {}", path.display()),
        ));
    }

    for name in ["raw_counts.tsv.gz", "raw_counts.tsv"] {
        let candidate = path.join(name);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    let mut candidates = Vec::new();
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let p = entry.path();
        if !p.is_file() {
            continue;
        }
        let Some(name) = p.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.ends_with("_raw_counts.tsv") || name.ends_with("_raw_counts.tsv.gz") {
            candidates.push(p);
        }
    }
    candidates.sort();
    candidates.into_iter().next().ok_or_else(|| {
        ScioError::new(
            ErrorCode::MissingFile,
            format!("expected raw_counts.tsv(.gz) in {}", path.display()),
        )
    })
}

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let resolved = resolve_bd_input_path(path)?;
    let mut md = crate::formats::dense::read_metadata(&resolved, strict)?;
    md.format = "bd_rhapsody_wta".to_string();
    Ok(md)
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    let resolved = resolve_bd_input_path(path)?;
    crate::formats::dense::read_matrix(&resolved, strict)
}
