use std::fmt::{Display, Formatter};

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

#[derive(Debug, Clone)]
pub struct ScioError {
    pub code: ErrorCode,
    pub message: String,
}

impl ScioError {
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

impl Display for ScioError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl std::error::Error for ScioError {}

impl From<std::io::Error> for ScioError {
    fn from(value: std::io::Error) -> Self {
        Self::new(ErrorCode::Io, value.to_string())
    }
}

pub type ScioResult<T> = Result<T, ScioError>;
