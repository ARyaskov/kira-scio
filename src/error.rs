use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCode {
    InvalidInputPath,
    MissingFile,
    UnsupportedFormat,
    ParseError,
    DimensionMismatch,
    ValidationError,
    FeatureDisabled,
    Io,
}

#[derive(Debug, Error)]
pub struct ScioError {
    pub code: ErrorCode,
    pub message: String,
    pub path: Option<PathBuf>,
    #[source]
    pub source: Option<Box<dyn std::error::Error + Send + Sync + 'static>>,
}

impl ScioError {
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            path: None,
            source: None,
        }
    }

    pub fn with_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.path = Some(path.into());
        self
    }

    pub fn with_source<E>(mut self, source: E) -> Self
    where
        E: std::error::Error + Send + Sync + 'static,
    {
        self.source = Some(Box::new(source));
        self
    }
}

impl Clone for ScioError {
    fn clone(&self) -> Self {
        Self {
            code: self.code,
            message: self.message.clone(),
            path: self.path.clone(),
            source: None,
        }
    }
}

impl std::fmt::Display for ScioError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.path.as_ref() {
            Some(p) => write!(f, "{:?} [{}]: {}", self.code, p.display(), self.message),
            None => write!(f, "{:?}: {}", self.code, self.message),
        }
    }
}

impl From<std::io::Error> for ScioError {
    fn from(value: std::io::Error) -> Self {
        let msg = value.to_string();
        ScioError::new(ErrorCode::Io, msg).with_source(value)
    }
}

pub type ScioResult<T> = Result<T, ScioError>;
