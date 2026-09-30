use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::detect::{DetectedFormat, detect_input_format};
use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{CanonicalData, InputMetadata, Marginals, MatrixStats, SoaCscMatrix};

/// Which feature modalities to keep when the source declares them.
///
/// 10x Feature Barcoding runs store antibody (`Antibody Capture`), CRISPR
/// and multiplexing features in the same matrix as genes; their counts are
/// on a different scale and distort gene-level statistics. Sources without a
/// feature type column are left untouched by every variant.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum FeatureTypeFilter {
    /// Keep every feature (default).
    #[default]
    All,
    /// Keep only features typed `Gene Expression` (case-insensitive).
    GeneExpression,
    /// Keep only features whose type matches one of these (case-insensitive).
    Only(Vec<String>),
}

impl FeatureTypeFilter {
    fn keeps(&self, feature_type: &str) -> bool {
        match self {
            Self::All => true,
            Self::GeneExpression => feature_type.eq_ignore_ascii_case("Gene Expression"),
            Self::Only(types) => types.iter().any(|t| t.eq_ignore_ascii_case(feature_type)),
        }
    }
}

/// Which matrix of an AnnData file to read.
///
/// Public h5ad files frequently carry a normalized `X` while the raw counts
/// live in `raw/X` or in a layer such as `layers/counts`. Check
/// `MatrixStats::is_integer` on the result if you need counts.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum H5adSource {
    /// The main matrix `/X` (default).
    #[default]
    X,
    /// `raw/X`, with gene labels taken from `raw/var`.
    RawX,
    /// `layers/<name>`, with gene labels from `/var`.
    Layer(String),
}

#[derive(Debug, Clone, Default)]
pub struct ReaderOptions {
    pub force_format: Option<DetectedFormat>,
    pub strict: bool,
    pub feature_types: FeatureTypeFilter,
    /// Only consulted for H5AD inputs.
    pub h5ad_source: H5adSource,
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
        Ok(self.load()?.0)
    }

    pub fn read_matrix(&self) -> ScioResult<SoaCscMatrix> {
        Ok(self.load()?.1)
    }

    pub fn read_all(&self) -> ScioResult<CanonicalData> {
        let (metadata, matrix) = self.load()?;
        Ok(CanonicalData { metadata, matrix })
    }

    /// Single parse of the input followed by the format-independent
    /// post-processing (feature-type filtering, marginals). Every public
    /// reader goes through here so labels, matrix and statistics always
    /// agree.
    fn load(&self) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
        let strict = self.options.strict;
        let (mut metadata, mut matrix) = match self.detected_format()? {
            DetectedFormat::Mtx10x => crate::formats::mtx10x::read_mtx(&self.input, strict)?,
            DetectedFormat::BdRhapsodyWta => {
                crate::formats::bd_rhapsody::read_all(&self.input, strict)?
            }
            DetectedFormat::DenseTsvCsv => {
                crate::formats::dense::parse_dense_full(&self.input, strict)?
            }
            DetectedFormat::H5ad => {
                crate::formats::h5ad::read_all(&self.input, strict, &self.options.h5ad_source)?
            }
            DetectedFormat::Loom => {
                // Stub backend; surfaces FeatureDisabled.
                let m = crate::formats::loom::read_metadata(&self.input, strict)?;
                let mx = crate::formats::loom::read_matrix(&self.input, strict)?;
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
        self.apply_feature_filter(&mut metadata, &mut matrix)?;
        metadata.marginals = Marginals::from_matrix(&matrix);
        Ok((metadata, matrix))
    }

    fn apply_feature_filter(
        &self,
        metadata: &mut InputMetadata,
        matrix: &mut SoaCscMatrix,
    ) -> ScioResult<()> {
        if self.options.feature_types == FeatureTypeFilter::All {
            return Ok(());
        }
        let Some(types) = metadata.feature_types.as_ref() else {
            return Ok(());
        };
        let keep: Vec<bool> = types
            .iter()
            .map(|t| self.options.feature_types.keeps(t))
            .collect();
        let excluded = keep.iter().filter(|k| !**k).count();
        if excluded == 0 {
            return Ok(());
        }
        *matrix = matrix
            .retain_genes(&keep)
            .map_err(|e| e.with_path(self.input.clone()))?;
        let filter = |v: &Vec<String>| -> Vec<String> {
            v.iter()
                .zip(&keep)
                .filter(|(_, k)| **k)
                .map(|(s, _)| s.clone())
                .collect()
        };
        metadata.gene_ids = filter(&metadata.gene_ids);
        metadata.gene_symbols = filter(&metadata.gene_symbols);
        metadata.feature_types = Some(filter(types));
        metadata.n_genes = matrix.n_genes;
        metadata.stats = MatrixStats::from_matrix(matrix);
        metadata.report.excluded_features = excluded;
        Ok(())
    }
}
