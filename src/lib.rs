#![forbid(unsafe_code)]

pub mod api;
pub mod cache_paths;
pub mod cli;
pub mod detect;
pub mod error;
pub mod formats;
pub mod model;
pub mod normalize;

pub use api::{FeatureTypeFilter, H5adSource, Reader, ReaderOptions};
pub use detect::{DetectedFormat, detect_input_format};
pub use error::{ErrorCode, ScioError, ScioResult};
pub use formats::bd_rhapsody::resolve_bd_input_path;
pub use formats::mtx10x::{MtxDatasetPaths, detect_prefix, discover};
pub use model::{
    CanonicalData, CountMismatch, IngestReport, InputMetadata, Marginals, MatrixKind, MatrixStats,
    Provenance, SoaCscMatrix, TOP_FEATURES,
};
