use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::{ErrorCode, ScioError, ScioResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DetectedFormat {
    Mtx10x,
    BdRhapsodyWta,
    DenseTsvCsv,
    H5ad,
    /// 10x Genomics Cell Ranger `.h5` matrix (`*_feature_bc_matrix.h5`).
    TenxH5,
    Loom,
}

pub fn detect_input_format(path: &Path) -> ScioResult<DetectedFormat> {
    if path.is_dir() {
        if crate::formats::mtx10x::contains_mtx_dataset(path)? {
            return Ok(DetectedFormat::Mtx10x);
        }
        if crate::formats::bd_rhapsody::resolve_bd_input_path(path).is_ok() {
            return Ok(DetectedFormat::BdRhapsodyWta);
        }
        return Err(ScioError::new(
            ErrorCode::UnsupportedFormat,
            format!(
                "failed to auto-detect format in directory {}",
                path.display()
            ),
        )
        .with_path(path.to_path_buf()));
    }

    if !path.exists() {
        return Err(ScioError::new(
            ErrorCode::InvalidInputPath,
            format!("input path does not exist: {}", path.display()),
        )
        .with_path(path.to_path_buf()));
    }

    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    if name.ends_with(".h5ad") || name.ends_with(".h5ad.gz") {
        return Ok(DetectedFormat::H5ad);
    }
    if name.ends_with(".loom") {
        return Ok(DetectedFormat::Loom);
    }
    if name.ends_with(".h5") || name.ends_with(".hdf5") {
        return Ok(DetectedFormat::TenxH5);
    }

    if filename_is_bd_counts(&name) {
        return Ok(DetectedFormat::BdRhapsodyWta);
    }

    if name.ends_with(".mtx") || name.ends_with(".mtx.gz") {
        return Ok(DetectedFormat::Mtx10x);
    }

    if name.ends_with(".tsv")
        || name.ends_with(".tsv.gz")
        || name.ends_with(".csv")
        || name.ends_with(".csv.gz")
    {
        return Ok(DetectedFormat::DenseTsvCsv);
    }

    Ok(sniff_content_or_default(path)?)
}

/// BD Rhapsody count tables by file name (lower-cased by the caller).
///
/// Covers the legacy `raw_counts.tsv` convention and the Sequence Analysis
/// Pipeline output `<sample>_{RSEC,DBEC}_MolsPerCell.csv[.gz]`.
fn filename_is_bd_counts(name: &str) -> bool {
    name == "raw_counts.tsv"
        || name == "raw_counts.tsv.gz"
        || name.contains("_raw_counts.tsv")
        || name.contains(".raw_counts.tsv")
        || name.contains("_molspercell.csv")
        || name.contains("_molspercell.tsv")
}

fn sniff_content_or_default(path: &Path) -> ScioResult<DetectedFormat> {
    match sniff_first_lines(path, 8) {
        Ok(lines) => Ok(classify_sniffed_lines(&lines)),
        Err(_) => Ok(DetectedFormat::DenseTsvCsv),
    }
}

fn classify_sniffed_lines(lines: &[String]) -> DetectedFormat {
    let mut header: Option<&str> = None;
    let mut data_lines: Vec<&str> = Vec::new();
    let mut saw_comment = false;

    for line in lines {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if t.starts_with('#') || t.starts_with('%') {
            saw_comment = true;
            continue;
        }
        if header.is_none() {
            header = Some(t);
        } else {
            data_lines.push(t);
        }
    }

    let Some(header) = header else {
        return DetectedFormat::DenseTsvCsv;
    };

    // MTX header: "rows cols nnz" on a single whitespace-separated line.
    let header_cols = header.split_whitespace().count();
    if header_cols == 3
        && header
            .split_whitespace()
            .all(|t| t.parse::<usize>().is_ok())
    {
        return DetectedFormat::Mtx10x;
    }

    let header_tabs = header.bytes().filter(|c| *c == b'\t').count();
    let header_commas = header.bytes().filter(|c| *c == b',').count();
    if header_tabs == 0 && header_commas == 0 {
        return DetectedFormat::DenseTsvCsv;
    }
    let delim = if header_commas > header_tabs {
        ','
    } else {
        '\t'
    };

    // BD Rhapsody signals: the Sequence Analysis Pipeline writes `#` comment
    // lines above the header and a `Cell_Index` first column. Value type is
    // deliberately not used: BD tables hold integer molecule counts, and
    // fractional values only mean the table was normalized upstream, which
    // says nothing about its origin.
    let first_cell = header
        .split(delim)
        .next()
        .map(|c| c.trim().to_ascii_lowercase())
        .unwrap_or_default();
    if saw_comment || first_cell == "cell_index" {
        return DetectedFormat::BdRhapsodyWta;
    }

    DetectedFormat::DenseTsvCsv
}

fn sniff_first_lines(path: &Path, max_lines: usize) -> ScioResult<Vec<String>> {
    let file = File::open(path).map_err(|e| {
        ScioError::new(ErrorCode::Io, e.to_string())
            .with_path(path.to_path_buf())
            .with_source(e)
    })?;
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    let reader: Box<dyn Read> = if name.ends_with(".gz") {
        Box::new(GzDecoder::new(BufReader::with_capacity(8 * 1024, file)))
    } else {
        Box::new(file)
    };
    let buffered = BufReader::with_capacity(8 * 1024, reader);
    let mut out = Vec::with_capacity(max_lines);
    for (i, line) in buffered.lines().take(max_lines).enumerate() {
        match line {
            Ok(l) if i == 0 => out.push(crate::normalize::strip_bom(&l).to_string()),
            Ok(l) => out.push(l),
            Err(_) => break,
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn sniffer_promotes_bd_on_leading_comment() {
        assert_eq!(
            classify_sniffed_lines(&lines(&["#meta", "cellA\tcellB", "GENE\t1\t2"])),
            DetectedFormat::BdRhapsodyWta
        );
    }

    #[test]
    fn sniffer_promotes_bd_on_cell_index_header() {
        assert_eq!(
            classify_sniffed_lines(&lines(&["Cell_Index,GENE_A,GENE_B", "1,5,0"])),
            DetectedFormat::BdRhapsodyWta
        );
    }

    #[test]
    fn sniffer_keeps_dense_for_integers_and_for_floats() {
        assert_eq!(
            classify_sniffed_lines(&lines(&["cellA\tcellB", "GENE\t1\t2"])),
            DetectedFormat::DenseTsvCsv
        );
        // Fractional values mean "normalized", not "BD Rhapsody".
        assert_eq!(
            classify_sniffed_lines(&lines(&["gene\tC1\tC2", "G1\t0.53\t1.2"])),
            DetectedFormat::DenseTsvCsv
        );
    }

    #[test]
    fn sniffer_recognizes_mtx_header() {
        assert_eq!(
            classify_sniffed_lines(&lines(&["%%MatrixMarket", "2 2 2", "1 1 1"])),
            DetectedFormat::Mtx10x
        );
    }
}
