//! BD Rhapsody WTA dense reader — thin shim over `formats::dense`.
//!
//! Accepts a count table file or a pipeline output directory. Inside a
//! directory the lookup order is: legacy `raw_counts.tsv[.gz]`, any
//! `*_raw_counts.tsv[.gz]`, then Sequence Analysis Pipeline tables
//! `*_DBEC_MolsPerCell.csv[.gz]` (recommended by BD for downstream analysis),
//! `*_RSEC_MolsPerCell.csv[.gz]`, and finally any other `*_MolsPerCell.csv[.gz]`.
//! Ties within a class resolve by lexicographic file name for determinism.

use std::path::{Path, PathBuf};

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, ShapeProbe, SoaCscMatrix};

/// Priority class of a BD count-table file name; lower sorts first.
fn bd_candidate_rank(lower_name: &str) -> Option<u8> {
    let is_table = |stem: &str| {
        lower_name.ends_with(&format!("{stem}.csv"))
            || lower_name.ends_with(&format!("{stem}.csv.gz"))
            || lower_name.ends_with(&format!("{stem}.tsv"))
            || lower_name.ends_with(&format!("{stem}.tsv.gz"))
    };
    if lower_name == "raw_counts.tsv" || lower_name == "raw_counts.tsv.gz" {
        Some(0)
    } else if is_table("_raw_counts") {
        Some(1)
    } else if is_table("_dbec_molspercell") {
        Some(2)
    } else if is_table("_rsec_molspercell") {
        Some(3)
    } else if is_table("_molspercell") {
        Some(4)
    } else {
        None
    }
}

pub fn resolve_bd_input_path(path: &Path) -> ScioResult<PathBuf> {
    if path.is_file() {
        return Ok(path.to_path_buf());
    }
    if !path.is_dir() {
        return Err(ScioError::new(
            ErrorCode::InvalidInputPath,
            format!("invalid input path: {}", path.display()),
        )
        .with_path(path.to_path_buf()));
    }

    let mut candidates: Vec<(u8, String, PathBuf)> = Vec::new();
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let p = entry.path();
        if !p.is_file() {
            continue;
        }
        let Some(name) = p.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        let lower = name.to_ascii_lowercase();
        if let Some(rank) = bd_candidate_rank(&lower) {
            candidates.push((rank, name.to_string(), p));
        }
    }
    candidates.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
    candidates
        .into_iter()
        .next()
        .map(|(_, _, p)| p)
        .ok_or_else(|| {
            ScioError::new(
                ErrorCode::MissingFile,
                format!(
                    "expected raw_counts.tsv(.gz) or *_MolsPerCell.csv(.gz) in {}",
                    path.display()
                ),
            )
            .with_path(path.to_path_buf())
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

pub(crate) fn read_shape(path: &Path) -> ScioResult<ShapeProbe> {
    crate::formats::dense::read_shape(&resolve_bd_input_path(path)?)
}

pub(crate) fn read_all(path: &Path, strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    let resolved = resolve_bd_input_path(path)?;
    let (mut md, mx) = crate::formats::dense::parse_dense_full(&resolved, strict)?;
    md.format = "bd_rhapsody_wta".to_string();
    Ok((md, mx))
}

#[cfg(test)]
mod tests {
    use super::bd_candidate_rank;

    #[test]
    fn candidate_rank_prefers_legacy_then_dbec_then_rsec() {
        assert_eq!(bd_candidate_rank("raw_counts.tsv"), Some(0));
        assert_eq!(bd_candidate_rank("s1_raw_counts.tsv.gz"), Some(1));
        assert_eq!(bd_candidate_rank("s1_dbec_molspercell.csv"), Some(2));
        assert_eq!(bd_candidate_rank("s1_rsec_molspercell.csv.gz"), Some(3));
        assert_eq!(bd_candidate_rank("s1_molspercell.csv"), Some(4));
        assert_eq!(bd_candidate_rank("s1_bioproduct_stats.csv"), None);
        assert_eq!(bd_candidate_rank("matrix.mtx"), None);
    }
}
