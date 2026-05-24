//! 10x Matrix Market reader (`matrix.mtx[.gz]` + features/genes + barcodes).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use rustc_hash::FxHashSet;
use tracing::warn;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, MatrixStats, SoaCscMatrix};
use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol};

#[derive(Debug, Clone)]
pub struct MtxDatasetPaths {
    pub input_dir: PathBuf,
    pub prefix: Option<String>,
    pub matrix: PathBuf,
    pub features: Option<PathBuf>,
    pub genes: Option<PathBuf>,
    pub barcodes: Option<PathBuf>,
}

pub const SHARED_CACHE_BASENAME: &str = "kira-organelle.bin";

pub fn contains_mtx_dataset(path: &Path) -> ScioResult<bool> {
    if !path.is_dir() {
        return Ok(false);
    }
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else { continue };
        if name.contains("matrix.mtx")
            || name.contains("features.tsv")
            || name.contains("barcodes.tsv")
        {
            return Ok(true);
        }
    }
    Ok(false)
}

pub fn read_metadata(path: &Path, strict: bool) -> ScioResult<InputMetadata> {
    let (metadata, _) = read_mtx(path, strict)?;
    Ok(metadata)
}

pub fn read_matrix(path: &Path, strict: bool) -> ScioResult<SoaCscMatrix> {
    let (_, matrix) = read_mtx(path, strict)?;
    Ok(matrix)
}

/// Single-pass entry point used by [`Reader::read_all`].
pub(crate) fn read_mtx(path: &Path, strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    let ds = discover(path)?;
    let matrix_reader = open_maybe_gz_existing(&ds.matrix)?;
    let parsed = parse_matrix_market(matrix_reader, &ds.matrix, strict)?;

    let (mut gene_ids, mut gene_symbols) = if let Some(features_path) = ds.features.as_ref() {
        parse_features(features_path, strict)?
    } else if let Some(genes_path) = ds.genes.as_ref() {
        parse_features(genes_path, strict)?
    } else {
        let synth: Vec<String> = (0..parsed.n_genes)
            .map(|i| normalize_gene_id("", None, i))
            .collect();
        (synth.clone(), synth)
    };

    fix_length(
        &mut gene_ids,
        parsed.n_genes,
        strict,
        "gene",
        &ds.matrix,
        |i| normalize_gene_id("", None, i),
    )?;
    fix_length(
        &mut gene_symbols,
        parsed.n_genes,
        strict,
        "gene_symbol",
        &ds.matrix,
        |i| normalize_gene_symbol("", None, i),
    )?;

    let mut barcodes = if let Some(path) = ds.barcodes.as_ref() {
        parse_barcodes(path)?
    } else {
        (0..parsed.n_cells).map(normalize_barcode_idx).collect()
    };
    fix_length(
        &mut barcodes,
        parsed.n_cells,
        strict,
        "barcode",
        &ds.matrix,
        normalize_barcode_idx,
    )?;

    let stats = parsed.stats.clone();
    let matrix = parsed.into_csc();
    matrix.validate()?;

    let metadata = InputMetadata {
        format: "mtx10x".to_string(),
        n_cells: matrix.n_cells,
        n_genes: matrix.n_genes,
        gene_ids,
        gene_symbols,
        barcodes,
        stats,
    };

    Ok((metadata, matrix))
}

fn normalize_barcode_idx(idx: usize) -> String {
    crate::normalize::synth_barcode(idx)
}

fn fix_length(
    out: &mut Vec<String>,
    expected: usize,
    strict: bool,
    label: &'static str,
    source: &Path,
    synth: impl Fn(usize) -> String,
) -> ScioResult<()> {
    if out.len() == expected {
        return Ok(());
    }
    if strict {
        return Err(ScioError::new(
            ErrorCode::DimensionMismatch,
            format!(
                "{} vector length {} != matrix {} count {}",
                label,
                out.len(),
                label,
                expected
            ),
        )
        .with_path(source.to_path_buf()));
    }
    out.resize_with(expected, || "".to_string());
    for (i, v) in out.iter_mut().enumerate() {
        if v.is_empty() {
            *v = synth(i);
        }
    }
    Ok(())
}

/// Intermediate triplet form; `into_csc()` sorts into CSC.
struct ParsedMtx {
    n_genes: usize,
    n_cells: usize,
    /// `(col, row, value)` — col first so stable sort yields CSC.
    triplets: Vec<(u32, u32, f32)>,
    stats: MatrixStats,
}

impl ParsedMtx {
    fn into_csc(mut self) -> SoaCscMatrix {
        self.triplets
            .sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));

        let n_cells = self.n_cells;
        let nnz = self.triplets.len();
        let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
        let mut row_idx: Vec<u32> = Vec::with_capacity(nnz);
        let mut values: Vec<f32> = Vec::with_capacity(nnz);

        col_ptr.push(0);
        let mut cur_col: u32 = 0;
        for (col, row, val) in self.triplets {
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

        SoaCscMatrix {
            n_cells,
            n_genes: self.n_genes,
            col_ptr,
            row_idx,
            values,
        }
    }
}

fn parse_matrix_market(
    reader: BufReader<Box<dyn Read>>,
    source: &Path,
    strict: bool,
) -> ScioResult<ParsedMtx> {
    let mut n_rows = None::<usize>;
    let mut n_cols = None::<usize>;
    let mut triplets: Vec<(u32, u32, f32)> = Vec::new();

    let mut nnz = 0usize;
    let mut total = 0f64;
    let mut min = f32::MAX;
    let mut max = f32::MIN;

    let mut saw_nonfinite = false;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }

        if n_rows.is_none() {
            let mut header = t.split_whitespace();
            let r = header
                .next()
                .ok_or_else(|| header_err(source, line_no))?
                .parse::<usize>()
                .map_err(|_| {
                    ScioError::new(ErrorCode::ParseError, "invalid n_rows")
                        .with_path(source.to_path_buf())
                })?;
            let c = header
                .next()
                .ok_or_else(|| header_err(source, line_no))?
                .parse::<usize>()
                .map_err(|_| {
                    ScioError::new(ErrorCode::ParseError, "invalid n_cols")
                        .with_path(source.to_path_buf())
                })?;
            let _ = header.next(); // nnz hint, ignored
            n_rows = Some(r);
            n_cols = Some(c);
            continue;
        }

        let mut parts = t.split_whitespace();
        let row_1 = parts
            .next()
            .and_then(|s| s.parse::<usize>().ok())
            .ok_or_else(|| {
                ScioError::new(
                    ErrorCode::ParseError,
                    format!("malformed coordinate at line {}", line_no + 1),
                )
                .with_path(source.to_path_buf())
            })?;
        let col_1 = parts
            .next()
            .and_then(|s| s.parse::<usize>().ok())
            .ok_or_else(|| {
                ScioError::new(
                    ErrorCode::ParseError,
                    format!("malformed coordinate at line {}", line_no + 1),
                )
                .with_path(source.to_path_buf())
            })?;
        let val = parts
            .next()
            .and_then(|s| s.parse::<f32>().ok())
            .ok_or_else(|| {
                ScioError::new(
                    ErrorCode::ParseError,
                    format!("malformed value at line {}", line_no + 1),
                )
                .with_path(source.to_path_buf())
            })?;

        if row_1 == 0 || col_1 == 0 {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "MTX is 1-based; found zero index",
            )
            .with_path(source.to_path_buf()));
        }

        let row = row_1 - 1;
        let col = col_1 - 1;
        if row >= n_rows.unwrap_or(0) || col >= n_cols.unwrap_or(0) {
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!("index out of range at line {}", line_no + 1),
                )
                .with_path(source.to_path_buf()));
            }
            continue;
        }

        if !val.is_finite() {
            saw_nonfinite = true;
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!(
                        "non-finite value at line {} (use strict=false to ignore)",
                        line_no + 1
                    ),
                )
                .with_path(source.to_path_buf()));
            }
            continue;
        }

        if val != 0.0 {
            triplets.push((col as u32, row as u32, val));
            nnz += 1;
            total += val as f64;
            min = min.min(val);
            max = max.max(val);
        }
    }

    if saw_nonfinite && !strict {
        warn!(
            path = %source.display(),
            "MTX contained non-finite values; dropped (strict=false)"
        );
    }

    let rows = n_rows.ok_or_else(|| {
        ScioError::new(ErrorCode::ParseError, "missing MTX header")
            .with_path(source.to_path_buf())
    })?;
    let cols = n_cols.ok_or_else(|| {
        ScioError::new(ErrorCode::ParseError, "missing MTX header")
            .with_path(source.to_path_buf())
    })?;

    let sparsity = if rows == 0 || cols == 0 {
        1.0
    } else {
        1.0 - (nnz as f64 / ((rows * cols) as f64))
    };
    let stats = MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min } else { 0.0 },
        max_count: if nnz > 0 { max } else { 0.0 },
        sparsity,
    };

    Ok(ParsedMtx {
        n_genes: rows,
        n_cells: cols,
        triplets,
        stats,
    })
}

fn header_err(source: &Path, line_no: usize) -> ScioError {
    ScioError::new(
        ErrorCode::ParseError,
        format!("invalid MTX header line {}", line_no + 1),
    )
    .with_path(source.to_path_buf())
}

/// Returns `(gene_ids, gene_symbols)`; single-column rows reuse the id.
fn parse_features(path: &Path, strict: bool) -> ScioResult<(Vec<String>, Vec<String>)> {
    let reader = open_maybe_gz_existing(path)?;
    let mut ids = Vec::<String>::new();
    let mut symbols = Vec::<String>::new();
    let mut single_column_row: Option<usize> = None;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim_end_matches(['\r', '\n']);
        if t.trim().is_empty() {
            continue;
        }
        let mut it = t.splitn(3, '\t');
        let id = it.next().unwrap_or("").trim();
        let sym = it.next().map(|s| s.trim());
        if sym.is_none() {
            single_column_row.get_or_insert(line_no);
        }
        ids.push(normalize_gene_id(id, sym, ids.len()));
        symbols.push(normalize_gene_symbol(id, sym, symbols.len()));
    }

    // Strict-mode: a single-column genes.tsv row is treated as malformed.
    if strict
        && path
            .file_name()
            .and_then(|f| f.to_str())
            .map(|name| name.starts_with("genes.tsv"))
            .unwrap_or(false)
        && let Some(line_no) = single_column_row
    {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!(
                "legacy genes.tsv must have <id>\\t<symbol> per row (line {} has 1 column)",
                line_no + 1
            ),
        )
        .with_path(path.to_path_buf()));
    }

    Ok((ids, symbols))
}

fn parse_barcodes(path: &Path) -> ScioResult<Vec<String>> {
    let reader = open_maybe_gz_existing(path)?;
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        out.push(normalize_barcode(&line, i));
    }
    Ok(out)
}

pub fn discover(path: &Path) -> ScioResult<MtxDatasetPaths> {
    let input_dir = if path.is_dir() {
        path.to_path_buf()
    } else {
        path.parent()
            .ok_or_else(|| ScioError::new(ErrorCode::InvalidInputPath, "input has no parent"))?
            .to_path_buf()
    };

    // Single read_dir + in-memory lookup, instead of ≤8 stat() probes per file.
    use rustc_hash::FxHashMap;
    let mut entries: FxHashMap<String, PathBuf> = FxHashMap::default();
    let mut prefixes: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for entry in std::fs::read_dir(&input_dir)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else { continue };
        entries.insert(name.to_string(), entry.path());
        if let Some(p) = extract_prefix(name) {
            prefixes.insert(p.to_string());
        }
    }
    if prefixes.len() > 1 {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!("multiple dataset prefixes in {}", input_dir.display()),
        )
        .with_path(input_dir));
    }
    let prefix = prefixes.into_iter().next();

    let lookup = |base: &str| -> Option<PathBuf> {
        let mut candidates: Vec<String> = Vec::with_capacity(6);
        if let Some(p) = prefix.as_deref() {
            candidates.push(format!("{p}_{base}"));
            candidates.push(format!("{p}_{base}.gz"));
            candidates.push(format!("{p}.{base}"));
            candidates.push(format!("{p}.{base}.gz"));
        }
        candidates.push(base.to_string());
        candidates.push(format!("{base}.gz"));
        for c in candidates {
            if let Some(p) = entries.get(&c) {
                return Some(p.clone());
            }
        }
        None
    };

    let matrix = lookup("matrix.mtx").ok_or_else(|| {
        ScioError::new(
            ErrorCode::MissingFile,
            format!("missing matrix.mtx(.gz) in {}", input_dir.display()),
        )
        .with_path(input_dir.clone())
    })?;

    let features = lookup("features.tsv");
    let genes = lookup("genes.tsv");
    let barcodes = lookup("barcodes.tsv");

    Ok(MtxDatasetPaths {
        input_dir,
        prefix,
        matrix,
        features,
        genes,
        barcodes,
    })
}

pub fn detect_prefix(input_dir: &Path) -> ScioResult<Option<String>> {
    let mut prefixes: BTreeSet<String> = BTreeSet::new();
    let mut seen = FxHashSet::default();
    for entry in std::fs::read_dir(input_dir)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else { continue };
        if !seen.insert(name.to_string()) {
            continue;
        }
        if let Some(p) = extract_prefix(name) {
            prefixes.insert(p.to_string());
        }
    }
    if prefixes.len() > 1 {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!("multiple dataset prefixes in {}", input_dir.display()),
        )
        .with_path(input_dir.to_path_buf()));
    }
    Ok(prefixes.into_iter().next())
}

pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    match prefix {
        Some(p) if !p.is_empty() => format!("{p}.{SHARED_CACHE_BASENAME}"),
        _ => SHARED_CACHE_BASENAME.to_string(),
    }
}

fn extract_prefix(name: &str) -> Option<&str> {
    // Longest-suffix-first so ".matrix.mtx.gz" beats ".matrix.mtx".
    const SUFFIXES: &[&str] = &[
        "_matrix.mtx.gz",
        ".matrix.mtx.gz",
        "_features.tsv.gz",
        ".features.tsv.gz",
        "_barcodes.tsv.gz",
        ".barcodes.tsv.gz",
        "_genes.tsv.gz",
        ".genes.tsv.gz",
        "_matrix.mtx",
        ".matrix.mtx",
        "_features.tsv",
        ".features.tsv",
        "_barcodes.tsv",
        ".barcodes.tsv",
        "_genes.tsv",
        ".genes.tsv",
    ];
    for suffix in SUFFIXES {
        if let Some(prefix) = name.strip_suffix(suffix)
            && !prefix.is_empty()
        {
            return Some(prefix);
        }
    }
    None
}

pub fn candidate_path(input_dir: &Path, prefix: Option<&str>, name: &str) -> PathBuf {
    match prefix {
        Some(p) if !p.is_empty() => {
            let underscore = input_dir.join(format!("{p}_{name}"));
            if exists_plain_or_gz(&underscore) {
                return underscore;
            }
            let dotted = input_dir.join(format!("{p}.{name}"));
            if exists_plain_or_gz(&dotted) {
                return dotted;
            }
            underscore
        }
        _ => input_dir.join(name),
    }
}

pub fn choose_existing(path: &Path) -> Option<PathBuf> {
    if path.exists() {
        return Some(path.to_path_buf());
    }
    let gz = gz_path(path);
    if gz.exists() {
        return Some(gz);
    }
    None
}

pub fn exists_plain_or_gz(path: &Path) -> bool {
    path.exists() || gz_path(path).exists()
}

pub fn gz_path(path: &Path) -> PathBuf {
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        path.with_extension(format!("{ext}.gz"))
    } else {
        path.with_extension("gz")
    }
}

pub fn open_maybe_gz_existing(path: &Path) -> ScioResult<BufReader<Box<dyn Read>>> {
    let existing = choose_existing(path).ok_or_else(|| {
        ScioError::new(
            ErrorCode::MissingFile,
            format!("missing file: {}", path.display()),
        )
        .with_path(path.to_path_buf())
    })?;
    let f = File::open(&existing).map_err(|e| {
        ScioError::new(ErrorCode::Io, e.to_string())
            .with_path(existing.clone())
            .with_source(e)
    })?;
    if existing
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        // Buffer the raw file so GzDecoder gets larger reads.
        let buffered_file = BufReader::with_capacity(64 * 1024, f);
        Ok(BufReader::new(Box::new(GzDecoder::new(buffered_file))))
    } else {
        Ok(BufReader::new(Box::new(f)))
    }
}
