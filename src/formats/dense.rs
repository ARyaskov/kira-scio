use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, MatrixStats, SoaCscMatrix};
use crate::normalize::{normalize_barcode, normalize_gene_symbol};

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let parsed = parse_dense(path, strict)?;
    Ok(InputMetadata {
        format: "dense_tsv_csv".to_string(),
        n_cells: parsed.barcodes.len(),
        n_genes: parsed.genes.len(),
        gene_ids: parsed.genes.clone(),
        gene_symbols: parsed.genes,
        barcodes: parsed.barcodes,
        stats: parsed.stats,
    })
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    Ok(parse_dense(path, strict)?.matrix)
}

#[derive(Debug)]
struct ParsedDense {
    genes: Vec<String>,
    barcodes: Vec<String>,
    matrix: SoaCscMatrix,
    stats: MatrixStats,
}

fn parse_dense(path: &Path, strict: bool) -> ScioResult<ParsedDense> {
    let mut genes = Vec::<String>::new();
    let mut barcodes = Vec::<String>::new();
    let mut columns = Vec::<Vec<(usize, f32)>>::new();

    let mut duplicate_genes = BTreeSet::new();
    let mut duplicate_barcodes = BTreeSet::new();
    let mut seen_genes = BTreeSet::new();
    let mut seen_barcodes = BTreeSet::new();

    let mut nnz = 0usize;
    let mut total = 0f64;
    let mut min_count = f32::MAX;
    let mut max_count = f32::MIN;

    let reader = open_maybe_gz(path)?;
    let mut header_parsed = false;
    let mut cell_major = false;
    let mut delim = '\t';

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.trim().is_empty() || trimmed.trim_start().starts_with('#') {
            continue;
        }

        if !header_parsed {
            delim = detect_delimiter(trimmed);
            let header = split_line(trimmed, delim);
            cell_major = first_col_looks_like_cell_header(header.first().map(|s| s.as_str()));

            if cell_major {
                genes = header
                    .iter()
                    .skip(1)
                    .enumerate()
                    .map(|(i, s)| normalize_gene_symbol(s, None, i))
                    .collect();
                if genes.is_empty() {
                    return Err(ScioError::new(
                        ErrorCode::ParseError,
                        "dense header has no genes",
                    ));
                }
            } else {
                let mut raw =
                    if first_col_looks_like_gene_header(header.first().map(|s| s.as_str())) {
                        header.into_iter().skip(1).collect::<Vec<_>>()
                    } else {
                        header
                    };
                if raw.is_empty() {
                    return Err(ScioError::new(
                        ErrorCode::ParseError,
                        "dense header has no barcodes",
                    ));
                }
                barcodes = raw
                    .drain(..)
                    .enumerate()
                    .map(|(i, s)| normalize_barcode(&s, i))
                    .collect();
                columns = vec![Vec::new(); barcodes.len()];
            }
            header_parsed = true;
            continue;
        }

        let mut parts = split_line(trimmed, delim);
        if cell_major {
            let expected = genes.len() + 1;
            if parts.len() + 1 == expected && !strict {
                parts.push(String::new());
            }
            if parts.len() != expected {
                return Err(ScioError::new(
                    ErrorCode::ParseError,
                    format!(
                        "line {}: expected {} columns, got {}",
                        line_no + 1,
                        expected,
                        parts.len()
                    ),
                ));
            }
            let barcode = normalize_barcode(&parts[0], barcodes.len());
            if !seen_barcodes.insert(barcode.clone()) {
                duplicate_barcodes.insert(barcode.clone());
            }
            barcodes.push(barcode);
            let mut col = Vec::new();
            for (g_idx, token) in parts.iter().skip(1).enumerate() {
                let value = parse_value(token, path, line_no + 1, g_idx + 2)?;
                if value == 0.0 {
                    continue;
                }
                nnz += 1;
                total += value as f64;
                min_count = min_count.min(value);
                max_count = max_count.max(value);
                col.push((g_idx, value));
            }
            columns.push(col);
        } else {
            let expected = barcodes.len() + 1;
            if parts.len() + 1 == expected && !strict {
                parts.push(String::new());
            }
            if parts.len() != expected {
                return Err(ScioError::new(
                    ErrorCode::ParseError,
                    format!(
                        "line {}: expected {} columns, got {}",
                        line_no + 1,
                        expected,
                        parts.len()
                    ),
                ));
            }
            let gene = normalize_gene_symbol(&parts[0], None, genes.len());
            if !seen_genes.insert(gene.clone()) {
                duplicate_genes.insert(gene.clone());
            }
            let row_idx = genes.len();
            genes.push(gene);
            for (c_idx, token) in parts.iter().skip(1).enumerate() {
                let value = parse_value(token, path, line_no + 1, c_idx + 2)?;
                if value == 0.0 {
                    continue;
                }
                nnz += 1;
                total += value as f64;
                min_count = min_count.min(value);
                max_count = max_count.max(value);
                columns[c_idx].push((row_idx, value));
            }
        }
    }

    if genes.is_empty() || barcodes.is_empty() {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!("empty dense matrix in {}", path.display()),
        ));
    }

    if !duplicate_genes.is_empty() && strict {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!(
                "duplicate gene symbols: {}",
                duplicate_genes.into_iter().collect::<Vec<_>>().join(",")
            ),
        ));
    }
    if !duplicate_barcodes.is_empty() && strict {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!(
                "duplicate barcodes: {}",
                duplicate_barcodes.into_iter().collect::<Vec<_>>().join(",")
            ),
        ));
    }

    let mut col_ptr = Vec::with_capacity(columns.len() + 1);
    let mut row_idx = Vec::with_capacity(nnz);
    let mut values = Vec::with_capacity(nnz);
    col_ptr.push(0);
    for col in &mut columns {
        col.sort_by_key(|(r, _)| *r);
        for (r, v) in col.iter() {
            row_idx.push(*r);
            values.push(*v);
        }
        col_ptr.push(row_idx.len());
    }

    let matrix = SoaCscMatrix {
        n_cells: barcodes.len(),
        n_genes: genes.len(),
        col_ptr,
        row_idx,
        values,
    };
    matrix.validate()?;

    let sparsity = 1.0 - (nnz as f64 / ((matrix.n_cells * matrix.n_genes) as f64));
    let stats = MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min_count } else { 0.0 },
        max_count: if nnz > 0 { max_count } else { 0.0 },
        sparsity,
    };

    Ok(ParsedDense {
        genes,
        barcodes,
        matrix,
        stats,
    })
}

fn parse_value(token: &str, path: &Path, line: usize, col: usize) -> ScioResult<f32> {
    let t = token.trim();
    if t.is_empty() {
        return Ok(0.0);
    }
    let mut val: f32 = t.parse().map_err(|_| {
        ScioError::new(
            ErrorCode::ParseError,
            format!(
                "{} line {} col {}: invalid numeric token `{}`",
                path.display(),
                line,
                col,
                t
            ),
        )
    })?;
    if !val.is_finite() {
        val = 0.0;
    }
    Ok(val)
}

fn open_maybe_gz(path: &Path) -> ScioResult<BufReader<Box<dyn Read>>> {
    if path
        .extension()
        .and_then(|v| v.to_str())
        .map(|v| v.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        let file = File::open(path)?;
        Ok(BufReader::new(Box::new(GzDecoder::new(file))))
    } else {
        let file = File::open(path)?;
        Ok(BufReader::new(Box::new(file)))
    }
}

fn detect_delimiter(line: &str) -> char {
    let commas = line.chars().filter(|c| *c == ',').count();
    let tabs = line.chars().filter(|c| *c == '\t').count();
    if commas > tabs { ',' } else { '\t' }
}

fn split_line(line: &str, delim: char) -> Vec<String> {
    line.split(delim).map(|s| s.trim().to_string()).collect()
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
