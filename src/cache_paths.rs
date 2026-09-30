//! Path helpers for the Kira shared single-cell cache (`kira-organelle.bin`)
//! that sits next to an MTX dataset.
//!
//! These are conventions of the Kira tool stack, not part of the reader
//! API; they live in their own module so the crate root only exposes
//! ingestion.

use std::path::{Path, PathBuf};

use crate::error::ScioResult;

/// Base name of the shared cache file written next to a dataset.
pub const SHARED_CACHE_BASENAME: &str = "kira-organelle.bin";

/// `<prefix>.kira-organelle.bin` for a prefixed dataset, else the base name.
pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    match prefix {
        Some(p) if !p.is_empty() => format!("{p}.{SHARED_CACHE_BASENAME}"),
        _ => SHARED_CACHE_BASENAME.to_string(),
    }
}

/// The single dataset prefix used by the MEX files in `input_dir`, if any.
/// Errors when the directory mixes several prefixes.
pub fn detect_prefix(input_dir: &Path) -> ScioResult<Option<String>> {
    crate::formats::mtx10x::detect_prefix(input_dir)
}

/// `<dir>/<prefix>_<name>` when it (or its `.gz`) exists, else
/// `<dir>/<prefix>.<name>` when that exists, else the underscore form;
/// `<dir>/<name>` without a prefix.
pub fn candidate_path(input_dir: &Path, prefix: Option<&str>, name: &str) -> PathBuf {
    match prefix {
        Some(p) if !p.is_empty() => {
            let underscore = input_dir.join(format!("{p}_{name}"));
            if exists_plain_or_gz(&underscore) {
                return underscore;
            }
            let dotted = input_dir.join(format!("{p}.{name}"));
            if exists_plain_or_gz(&dotted) {
                return dotted;
            }
            underscore
        }
        _ => input_dir.join(name),
    }
}

/// `path` if it exists, else `path.gz` if that exists.
pub fn choose_existing(path: &Path) -> Option<PathBuf> {
    if path.exists() {
        return Some(path.to_path_buf());
    }
    let gz = gz_path(path);
    if gz.exists() {
        return Some(gz);
    }
    None
}

pub fn exists_plain_or_gz(path: &Path) -> bool {
    path.exists() || gz_path(path).exists()
}

/// Appends `.gz` to the file name, keeping any existing extension.
pub fn gz_path(path: &Path) -> PathBuf {
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        path.with_extension(format!("{ext}.gz"))
    } else {
        path.with_extension("gz")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_filename_follows_prefix() {
        assert_eq!(resolve_shared_cache_filename(None), "kira-organelle.bin");
        assert_eq!(
            resolve_shared_cache_filename(Some("")),
            "kira-organelle.bin"
        );
        assert_eq!(
            resolve_shared_cache_filename(Some("S1")),
            "S1.kira-organelle.bin"
        );
    }

    #[test]
    fn gz_path_keeps_the_original_extension() {
        assert_eq!(
            gz_path(Path::new("/d/matrix.mtx")),
            PathBuf::from("/d/matrix.mtx.gz")
        );
        assert_eq!(gz_path(Path::new("/d/noext")), PathBuf::from("/d/noext.gz"));
    }

    #[test]
    fn candidate_and_existing_resolution() {
        let dir = std::env::temp_dir().join(format!(
            "kira_scio_cache_paths_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("S1.matrix.mtx.gz"), b"").unwrap();

        // Dotted form exists (as .gz), underscore form does not.
        assert_eq!(
            candidate_path(&dir, Some("S1"), "matrix.mtx"),
            dir.join("S1.matrix.mtx")
        );
        assert_eq!(
            choose_existing(&dir.join("S1.matrix.mtx")),
            Some(dir.join("S1.matrix.mtx.gz"))
        );
        assert_eq!(choose_existing(&dir.join("absent.mtx")), None);
        // No prefix: plain join, whether or not it exists.
        assert_eq!(
            candidate_path(&dir, None, "matrix.mtx"),
            dir.join("matrix.mtx")
        );
        // Unknown prefix falls back to the underscore form.
        assert_eq!(
            candidate_path(&dir, Some("Z"), "matrix.mtx"),
            dir.join("Z_matrix.mtx")
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
