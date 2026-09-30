use std::path::Path;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, ShapeProbe, SoaCscMatrix};

pub fn read_metadata(path: &Path, _strict: bool) -> ScioResult<InputMetadata> {
    let _ = path;
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "loom reader is optional and not yet enabled in this build",
    ))
}

pub(crate) fn read_shape(path: &Path) -> ScioResult<ShapeProbe> {
    let _ = path;
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "loom reader is optional and not yet enabled in this build",
    ))
}

pub fn read_matrix(path: &Path, _strict: bool) -> ScioResult<SoaCscMatrix> {
    let _ = path;
    Err(ScioError::new(
        ErrorCode::FeatureDisabled,
        "loom reader is optional and not yet enabled in this build",
    ))
}
