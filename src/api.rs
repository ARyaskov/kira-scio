use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::detect::{DetectedFormat, detect_input_format};
use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{CanonicalData, InputMetadata, SoaCscMatrix};

#[derive(Debug, Clone, Default)]
pub struct ReaderOptions {
    pub force_format: Option<DetectedFormat>,
    pub strict: bool,
}

#[derive(Debug)]
pub struct Reader {
    input: PathBuf,
    options: ReaderOptions,
    detected: OnceLock<DetectedFormat>,
}

impl Clone for Reader {
    fn clone(&self) -> Self {
        Self {
            input: self.input.clone(),
            options: self.options.clone(),
            detected: OnceLock::new(),
        }
    }
}

impl Reader {
    pub fn new(input: impl AsRef<Path>) -> Self {
        Self {
            input: input.as_ref().to_path_buf(),
            options: ReaderOptions {
                strict: true,
                ..ReaderOptions::default()
            },
            detected: OnceLock::new(),
        }
    }

    pub fn with_options(input: impl AsRef<Path>, options: ReaderOptions) -> Self {
        Self {
            input: input.as_ref().to_path_buf(),
            options,
            detected: OnceLock::new(),
        }
    }

    pub fn detected_format(&self) -> ScioResult<DetectedFormat> {
        if let Some(fmt) = self.options.force_format {
            return Ok(fmt);
        }
        if let Some(fmt) = self.detected.get() {
            return Ok(*fmt);
        }
        let fmt = detect_input_format(&self.input)?;
        let _ = self.detected.set(fmt);
        Ok(fmt)
    }

    pub fn read_metadata(&self) -> ScioResult<InputMetadata> {
        match self.detected_format()? {
            DetectedFormat::Mtx10x => {
                crate::formats::mtx10x::read_metadata(&self.input, self.options.strict)
            }
            DetectedFormat::BdRhapsodyWta => {
                crate::formats::bd_rhapsody::read_metadata(&self.input, self.options.strict)
            }
            DetectedFormat::DenseTsvCsv => {
                crate::formats::dense::read_metadata(&self.input, self.options.strict)
            }
            DetectedFormat::H5ad => {
                crate::formats::h5ad::read_metadata(&self.input, self.options.strict)
            }
            DetectedFormat::Loom => {
                crate::formats::loom::read_metadata(&self.input, self.options.strict)
            }
        }
    }

    pub fn read_matrix(&self) -> ScioResult<SoaCscMatrix> {
        match self.detected_format()? {
            DetectedFormat::Mtx10x => {
                crate::formats::mtx10x::read_matrix(&self.input, self.options.strict)
            }
            DetectedFormat::BdRhapsodyWta => {
                crate::formats::bd_rhapsody::read_matrix(&self.input, self.options.strict)
            }
            DetectedFormat::DenseTsvCsv => {
                crate::formats::dense::read_matrix(&self.input, self.options.strict)
            }
            DetectedFormat::H5ad => {
                crate::formats::h5ad::read_matrix(&self.input, self.options.strict)
            }
            DetectedFormat::Loom => {
                crate::formats::loom::read_matrix(&self.input, self.options.strict)
            }
        }
    }

    pub fn read_all(&self) -> ScioResult<CanonicalData> {
        let (metadata, matrix) = match self.detected_format()? {
            DetectedFormat::Mtx10x => {
                crate::formats::mtx10x::read_mtx(&self.input, self.options.strict)?
            }
            DetectedFormat::BdRhapsodyWta => {
                crate::formats::bd_rhapsody::read_all(&self.input, self.options.strict)?
            }
            DetectedFormat::DenseTsvCsv => {
                crate::formats::dense::parse_dense_full(&self.input, self.options.strict)?
            }
            DetectedFormat::H5ad => {
                crate::formats::h5ad::read_all(&self.input, self.options.strict)?
            }
            DetectedFormat::Loom => {
                // Stub backend; surfaces FeatureDisabled.
                let m = crate::formats::loom::read_metadata(&self.input, self.options.strict)?;
                let mx = crate::formats::loom::read_matrix(&self.input, self.options.strict)?;
                (m, mx)
            }
        };
        if metadata.n_cells != matrix.n_cells || metadata.n_genes != matrix.n_genes {
            return Err(ScioError::new(
                ErrorCode::DimensionMismatch,
                format!(
                    "metadata/matrix dimensions mismatch: metadata={}x{}, matrix={}x{}",
                    metadata.n_cells, metadata.n_genes, matrix.n_cells, matrix.n_genes
                ),
            )
            .with_path(self.input.clone()));
        }
        Ok(CanonicalData { metadata, matrix })
    }
}
