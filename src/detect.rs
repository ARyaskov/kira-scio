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

    if filename_is_bd_raw_counts(&name) {
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

fn filename_is_bd_raw_counts(name: &str) -> bool {
    name == "raw_counts.tsv"
        || name == "raw_counts.tsv.gz"
        || name.ends_with("_raw_counts.tsv")
        || name.ends_with("_raw_counts.tsv.gz")
        || name.contains("_raw_counts.tsv")
        || name.contains(".raw_counts.tsv")
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
    let delim: u8 = if header_commas > header_tabs {
        b','
    } else {
        b'\t'
    };

    // Heuristic 1: a leading comment line is the strongest BD Rhapsody signal.
    if saw_comment {
        return DetectedFormat::BdRhapsodyWta;
    }

    // Heuristic 2: at least one data row with float values (e.g. ".5", "1.0")
    // suggests normalized BD Rhapsody output. Plain integers default to
    // generic dense TSV so 10x/MEX-shaped TSVs don't get mis-promoted.
    let any_float = data_lines.iter().any(|line| has_float_value(line, delim));
    if any_float {
        return DetectedFormat::BdRhapsodyWta;
    }

    DetectedFormat::DenseTsvCsv
}

fn has_float_value(line: &str, delim: u8) -> bool {
    let mut first = true;
    for token in line.split(delim as char) {
        if first {
            // Skip the gene/cell label column.
            first = false;
            continue;
        }
        let t = token.trim();
        if t.is_empty() {
            continue;
        }
        if t.contains('.') || t.contains('e') || t.contains('E') {
            if t.parse::<f64>().is_ok() {
                return true;
            }
        }
    }
    false
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
    for line in buffered.lines().take(max_lines) {
        match line {
            Ok(l) => out.push(l),
            Err(_) => break,
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sniffer_promotes_bd_on_leading_comment() {
        let lines = vec![
            "#meta".to_string(),
            "cellA\tcellB".to_string(),
            "GENE\t1.0\t2.5".to_string(),
        ];
        assert_eq!(
            classify_sniffed_lines(&lines),
            DetectedFormat::BdRhapsodyWta
        );
    }

    #[test]
    fn sniffer_defaults_to_dense_for_integers() {
        let lines = vec!["cellA\tcellB".to_string(), "GENE\t1\t2".to_string()];
        assert_eq!(classify_sniffed_lines(&lines), DetectedFormat::DenseTsvCsv);
    }

    #[test]
    fn sniffer_recognizes_mtx_header() {
        let lines = vec![
            "%%MatrixMarket".to_string(),
            "2 2 2".to_string(),
            "1 1 1".to_string(),
        ];
        assert_eq!(classify_sniffed_lines(&lines), DetectedFormat::Mtx10x);
    }
}
