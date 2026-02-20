use std::path::Path;

use crate::api::Reader;
use crate::error::ScioResult;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InspectCommand {
    Inspect,
    Stats,
    Validate,
}

pub fn run(command: InspectCommand, input: &Path) -> ScioResult<String> {
    let reader = Reader::new(input);
    match command {
        InspectCommand::Inspect => {
            let md = reader.read_metadata()?;
            Ok(format!(
                "format={} n_cells={} n_genes={} nnz={}",
                md.format, md.n_cells, md.n_genes, md.stats.nnz
            ))
        }
        InspectCommand::Stats => {
            let md = reader.read_metadata()?;
            Ok(format!(
                "sparsity={:.6} total_counts={:.3} min={:.3} max={:.3}",
                md.stats.sparsity, md.stats.total_counts, md.stats.min_count, md.stats.max_count
            ))
        }
        InspectCommand::Validate => {
            let data = reader.read_all()?;
            data.matrix.validate()?;
            Ok("ok".to_string())
        }
    }
}
