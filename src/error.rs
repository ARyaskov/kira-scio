use std::path::PathBuf;

use thiserror::Error;

/// Stable, machine-readable failure category. New variants may be added in
/// a minor release; match with a wildcard arm. [`ErrorCode::as_str`] gives
/// the snake_case code that stays fixed across versions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
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

impl ErrorCode {
    /// Stable snake_case code for logs, exit statuses and JSON.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InvalidInputPath => "invalid_input_path",
            Self::MissingFile => "missing_file",
            Self::UnsupportedFormat => "unsupported_format",
            Self::ParseError => "parse_error",
            Self::DimensionMismatch => "dimension_mismatch",
            Self::ValidationError => "validation_error",
            Self::FeatureDisabled => "feature_disabled",
            Self::Io => "io",
        }
    }

    /// All codes, for exhaustive documentation and tests.
    pub const ALL: &'static [ErrorCode] = &[
        Self::InvalidInputPath,
        Self::MissingFile,
        Self::UnsupportedFormat,
        Self::ParseError,
        Self::DimensionMismatch,
        Self::ValidationError,
        Self::FeatureDisabled,
        Self::Io,
    ];
}

impl std::fmt::Display for ErrorCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
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
            Some(p) => write!(f, "{} [{}]: {}", self.code, p.display(), self.message),
            None => write!(f, "{}: {}", self.code, self.message),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable_snake_case_and_unique() {
        let mut seen = std::collections::HashSet::new();
        for code in ErrorCode::ALL {
            let s = code.as_str();
            assert!(
                s.bytes().all(|b| b.is_ascii_lowercase() || b == b'_'),
                "{s}"
            );
            assert!(seen.insert(s), "duplicate code {s}");
            assert_eq!(code.to_string(), s);
        }
        assert_eq!(ErrorCode::ParseError.as_str(), "parse_error");
    }

    #[test]
    fn display_carries_code_path_and_message() {
        let e = ScioError::new(ErrorCode::MissingFile, "no matrix").with_path("/d/x");
        assert_eq!(e.to_string(), "missing_file [/d/x]: no matrix");
        let e = ScioError::new(ErrorCode::Io, "boom");
        assert_eq!(e.to_string(), "io: boom");
        let cloned = e.clone();
        assert_eq!(cloned.code, ErrorCode::Io);
        assert!(cloned.source.is_none());
    }
}
