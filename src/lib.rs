#![forbid(unsafe_code)]

pub mod api;
pub mod cli;
pub mod detect;
pub mod error;
pub mod formats;
pub mod model;
pub mod normalize;

pub use api::{Reader, ReaderOptions};
pub use detect::{DetectedFormat, detect_input_format};
pub use error::{ErrorCode, ScioError, ScioResult};
pub use formats::bd_rhapsody::resolve_bd_input_path;
pub use formats::mtx10x::{
    MtxDatasetPaths, SHARED_CACHE_BASENAME, candidate_path, choose_existing, detect_prefix,
    discover, exists_plain_or_gz, gz_path, open_maybe_gz_existing, resolve_shared_cache_filename,
};
pub use model::{CanonicalData, InputMetadata, MatrixStats, SoaCscMatrix};
