//! Dense TSV/CSV reader (BD Rhapsody WTA and similar).
//!
//! Auto-detects row-major vs cell-major orientation. Duplicates are kept in
//! insertion order; a `tracing::warn!` is emitted regardless of strict mode.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;
use rustc_hash::FxHashSet;
use tracing::warn;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, MatrixStats, SoaCscMatrix};
use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol};

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let parsed = parse_dense(path, strict)?;
    Ok(InputMetadata {
        format: "dense_tsv_csv".to_string(),
        n_cells: parsed.barcodes.len(),
        n_genes: parsed.gene_ids.len(),
        gene_ids: parsed.gene_ids,
        gene_symbols: parsed.gene_symbols,
        barcodes: parsed.barcodes,
        stats: parsed.stats,
    })
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    Ok(parse_dense(path, strict)?.matrix)
}

/// Single-pass entry point used by [`Reader::read_all`].
pub(crate) fn parse_dense_full(
    path: &Path,
    strict: bool,
) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    let parsed = parse_dense(path, strict)?;
    let metadata = InputMetadata {
        format: "dense_tsv_csv".to_string(),
        n_cells: parsed.barcodes.len(),
        n_genes: parsed.gene_ids.len(),
        gene_ids: parsed.gene_ids,
        gene_symbols: parsed.gene_symbols,
        barcodes: parsed.barcodes,
        stats: parsed.stats,
    };
    Ok((metadata, parsed.matrix))
}

#[derive(Debug)]
struct ParsedDense {
    gene_ids: Vec<String>,
    gene_symbols: Vec<String>,
    barcodes: Vec<String>,
    matrix: SoaCscMatrix,
    stats: MatrixStats,
}

fn parse_dense(path: &Path, strict: bool) -> ScioResult<ParsedDense> {
    let mut gene_ids = Vec::<String>::new();
    let mut gene_symbols = Vec::<String>::new();
    let mut barcodes = Vec::<String>::new();
    // (col, row, value) triplets, sorted into CSC at the end.
    let mut triplets: Vec<(u32, u32, f32)> = Vec::new();

    let mut seen_genes = FxHashSet::default();
    let mut seen_barcodes = FxHashSet::default();
    let mut duplicate_genes: Vec<String> = Vec::new();
    let mut duplicate_barcodes: Vec<String> = Vec::new();

    let mut nnz = 0usize;
    let mut total = 0f64;
    let mut min_count = f32::MAX;
    let mut max_count = f32::MIN;
    let mut saw_nonfinite = false;

    let reader = open_maybe_gz(path)?;
    let mut header_parsed = false;
    let mut cell_major = false;
    let mut delim = '\t';

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.trim().is_empty() || trimmed.trim_start().starts_with('#') {
            continue;
        }

        if !header_parsed {
            delim = detect_delimiter(trimmed);
            let header = split_line(trimmed, delim);
            cell_major = first_col_looks_like_cell_header(header.first().map(|s| s.as_str()));

            if cell_major {
                gene_ids.reserve(header.len().saturating_sub(1));
                gene_symbols.reserve(header.len().saturating_sub(1));
                for (i, s) in header.iter().enumerate().skip(1) {
                    let id = normalize_gene_id(s, None, gene_ids.len());
                    let sym = normalize_gene_symbol(s, None, gene_symbols.len());
                    if !seen_genes.insert(id.clone()) {
                        duplicate_genes.push(id.clone());
                    }
                    gene_ids.push(id);
                    gene_symbols.push(sym);
                    let _ = i;
                }
                if gene_ids.is_empty() {
                    return Err(ScioError::new(
                        ErrorCode::ParseError,
                        "dense header has no genes",
                    )
                    .with_path(path.to_path_buf()));
                }
            } else {
                // Row-major. First cell: label ("gene"/"feature"), empty
                // (BD Rhapsody `\tC1\tC2`), or already a barcode.
                let starts_with_label =
                    first_col_looks_like_gene_header(header.first().map(|s| s.as_str()));
                let starts_with_empty = header
                    .first()
                    .map(|s| s.is_empty())
                    .unwrap_or(false);
                let mut raw = if starts_with_label || starts_with_empty {
                    header.into_iter().skip(1).collect::<Vec<_>>()
                } else {
                    header
                };
                if raw.is_empty() {
                    return Err(ScioError::new(
                        ErrorCode::ParseError,
                        "dense header has no barcodes",
                    )
                    .with_path(path.to_path_buf()));
                }
                barcodes.reserve(raw.len());
                for (i, s) in raw.drain(..).enumerate() {
                    let b = normalize_barcode(&s, i);
                    if !seen_barcodes.insert(b.clone()) {
                        duplicate_barcodes.push(b.clone());
                    }
                    barcodes.push(b);
                }
            }
            header_parsed = true;
            continue;
        }

        if cell_major {
            let expected = gene_ids.len() + 1;
            let raw_parts = count_fields(trimmed, delim);
            // Tolerate a missing trailing column (common BD writer quirk).
            if raw_parts + 1 != expected && raw_parts != expected {
                return Err(ScioError::new(
                    ErrorCode::ParseError,
                    format!(
                        "line {}: expected {} columns, got {}",
                        line_no + 1,
                        expected,
                        raw_parts
                    ),
                )
                .with_path(path.to_path_buf()));
            }
            let mut tokens = trimmed.split(delim);
            let barcode_raw = tokens.next().unwrap_or("");
            let barcode = normalize_barcode(barcode_raw.trim(), barcodes.len());
            if !seen_barcodes.insert(barcode.clone()) {
                duplicate_barcodes.push(barcode.clone());
            }
            let col_idx = barcodes.len() as u32;
            barcodes.push(barcode);

            for g_idx in 0..gene_ids.len() {
                let token = tokens.next().unwrap_or("");
                let value = parse_value(token.trim(), path, line_no + 1, g_idx + 2, strict)?;
                match value {
                    ValueOutcome::Zero => {}
                    ValueOutcome::NonFinite => saw_nonfinite = true,
                    ValueOutcome::Finite(v) => {
                        triplets.push((col_idx, g_idx as u32, v));
                        nnz += 1;
                        total += v as f64;
                        min_count = min_count.min(v);
                        max_count = max_count.max(v);
                    }
                }
            }
        } else {
            let expected = barcodes.len() + 1;
            let raw_parts = count_fields(trimmed, delim);
            if raw_parts + 1 == expected || raw_parts == expected {
                // ok
            } else {
                return Err(ScioError::new(
                    ErrorCode::ParseError,
                    format!(
                        "line {}: expected {} columns, got {}",
                        line_no + 1,
                        expected,
                        raw_parts
                    ),
                )
                .with_path(path.to_path_buf()));
            }

            let mut tokens = trimmed.split(delim);
            let gene_raw = tokens.next().unwrap_or("");
            let gene_id = normalize_gene_id(gene_raw.trim(), None, gene_ids.len());
            let gene_symbol = normalize_gene_symbol(gene_raw.trim(), None, gene_symbols.len());
            if !seen_genes.insert(gene_id.clone()) {
                duplicate_genes.push(gene_id.clone());
            }
            let row_idx_value = gene_ids.len() as u32;
            gene_ids.push(gene_id);
            gene_symbols.push(gene_symbol);

            for c_idx in 0..barcodes.len() {
                let token = tokens.next().unwrap_or("");
                let value = parse_value(token.trim(), path, line_no + 1, c_idx + 2, strict)?;
                match value {
                    ValueOutcome::Zero => {}
                    ValueOutcome::NonFinite => saw_nonfinite = true,
                    ValueOutcome::Finite(v) => {
                        triplets.push((c_idx as u32, row_idx_value, v));
                        nnz += 1;
                        total += v as f64;
                        min_count = min_count.min(v);
                        max_count = max_count.max(v);
                    }
                }
            }
        }
    }

    if gene_ids.is_empty() || barcodes.is_empty() {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!("empty dense matrix in {}", path.display()),
        )
        .with_path(path.to_path_buf()));
    }

    if !duplicate_genes.is_empty() {
        warn!(path = %path.display(), duplicates = ?duplicate_genes, "duplicate gene symbols");
    }
    if !duplicate_barcodes.is_empty() {
        warn!(path = %path.display(), duplicates = ?duplicate_barcodes, "duplicate barcodes");
    }
    if saw_nonfinite && !strict {
        warn!(path = %path.display(), "non-finite values treated as zero");
    }

    triplets.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    let n_cells = barcodes.len();
    let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
    let mut row_idx: Vec<u32> = Vec::with_capacity(triplets.len());
    let mut values: Vec<f32> = Vec::with_capacity(triplets.len());
    col_ptr.push(0);
    let mut cur_col: u32 = 0;
    for (col, row, val) in triplets {
        while cur_col < col {
            col_ptr.push(row_idx.len() as u64);
            cur_col += 1;
        }
        row_idx.push(row);
        values.push(val);
    }
    while col_ptr.len() <= n_cells {
        col_ptr.push(row_idx.len() as u64);
    }

    let matrix = SoaCscMatrix {
        n_cells,
        n_genes: gene_ids.len(),
        col_ptr,
        row_idx,
        values,
    };
    matrix.validate()?;

    let total_cells = matrix.n_cells;
    let total_genes = matrix.n_genes;
    let sparsity = if total_cells == 0 || total_genes == 0 {
        1.0
    } else {
        1.0 - (nnz as f64 / ((total_cells * total_genes) as f64))
    };
    let stats = MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min_count } else { 0.0 },
        max_count: if nnz > 0 { max_count } else { 0.0 },
        sparsity,
    };

    Ok(ParsedDense {
        gene_ids,
        gene_symbols,
        barcodes,
        matrix,
        stats,
    })
}

enum ValueOutcome {
    Zero,
    NonFinite,
    Finite(f32),
}

fn parse_value(
    token: &str,
    path: &Path,
    line: usize,
    col: usize,
    strict: bool,
) -> ScioResult<ValueOutcome> {
    if token.is_empty() {
        return Ok(ValueOutcome::Zero);
    }
    let val: f32 = token.parse().map_err(|_| {
        ScioError::new(
            ErrorCode::ParseError,
            format!(
                "line {} col {}: invalid numeric token `{}`",
                line, col, token
            ),
        )
        .with_path(path.to_path_buf())
    })?;
    if !val.is_finite() {
        if strict {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!(
                    "line {} col {}: non-finite value `{}` (use strict=false to ignore)",
                    line, col, token
                ),
            )
            .with_path(path.to_path_buf()));
        }
        return Ok(ValueOutcome::NonFinite);
    }
    if val == 0.0 {
        Ok(ValueOutcome::Zero)
    } else {
        Ok(ValueOutcome::Finite(val))
    }
}

fn open_maybe_gz(path: &Path) -> ScioResult<BufReader<Box<dyn Read>>> {
    let file = File::open(path).map_err(|e| {
        ScioError::new(ErrorCode::Io, e.to_string())
            .with_path(path.to_path_buf())
            .with_source(e)
    })?;
    if path
        .extension()
        .and_then(|v| v.to_str())
        .map(|v| v.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        let buffered = BufReader::with_capacity(64 * 1024, file);
        Ok(BufReader::new(Box::new(GzDecoder::new(buffered))))
    } else {
        Ok(BufReader::new(Box::new(file)))
    }
}

fn detect_delimiter(line: &str) -> char {
    let commas = line.bytes().filter(|c| *c == b',').count();
    let tabs = line.bytes().filter(|c| *c == b'\t').count();
    if commas > tabs { ',' } else { '\t' }
}

// NOTE: quoted CSV (`"foo,bar",1,2`) is not handled; in-scope formats never quote.
fn split_line(line: &str, delim: char) -> Vec<String> {
    line.split(delim).map(|s| s.trim().to_string()).collect()
}

fn count_fields(line: &str, delim: char) -> usize {
    if line.is_empty() {
        0
    } else {
        line.split(delim).count()
    }
}

fn first_col_looks_like_gene_header(v: Option<&str>) -> bool {
    let Some(v) = v else { return false };
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "gene" | "genes" | "gene_symbol" | "genesymbol" | "symbol" | "feature" | "features"
    )
}

fn first_col_looks_like_cell_header(v: Option<&str>) -> bool {
    let Some(v) = v else { return false };
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "barcode" | "barcodes" | "cell" | "cell_id" | "cellid"
    )
}
