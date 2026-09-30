//! Shared helpers for integration tests.
//!
//! `TestDir` is a scratch directory that is deleted when dropped, so tests
//! leave nothing behind in the system temp dir even when they fail.
#![allow(dead_code)]

use std::io::Write;
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static COUNTER: AtomicU64 = AtomicU64::new(0);

/// A temporary directory removed on drop. Derefs to `Path`, so `d.join(..)`
/// and `Reader::new(&d)` work directly.
pub struct TestDir(PathBuf);

impl TestDir {
    pub fn path(&self) -> &Path {
        &self.0
    }

    /// Writes `content` to `<dir>/<name>` and returns the path.
    pub fn file(&self, name: &str, content: &str) -> PathBuf {
        let p = self.0.join(name);
        write(&p, content);
        p
    }

    /// Creates and returns a nested directory.
    pub fn subdir(&self, name: &str) -> PathBuf {
        let p = self.0.join(name);
        std::fs::create_dir_all(&p).unwrap();
        p
    }
}

impl Deref for TestDir {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for TestDir {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A fresh, uniquely named scratch directory under the system temp dir.
pub fn temp_dir(label: &str) -> TestDir {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "kira_scio_test_{}_{}_{}",
        std::process::id(),
        n,
        label
    ));
    std::fs::create_dir_all(&dir).unwrap();
    TestDir(dir)
}

pub fn write(path: &Path, content: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}
