use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{InputMetadata, MatrixStats, SoaCscMatrix};
use crate::normalize::{normalize_barcode, normalize_gene_symbol};

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

fn read_mtx(path: &Path, strict: bool) -> ScioResult<(InputMetadata, SoaCscMatrix)> {
    let ds = discover(path)?;
    let mut matrix_reader = open_maybe_gz_existing(&ds.matrix)?;

    let (n_genes, n_cells, mut columns, stats) =
        parse_matrix_market(&mut matrix_reader, &ds.matrix, strict)?;

    let mut gene_ids = if let Some(features_path) = ds.features.as_ref() {
        parse_features(features_path)?
    } else if let Some(genes_path) = ds.genes.as_ref() {
        parse_legacy_genes(genes_path)?
    } else {
        (0..n_genes)
            .map(|i| normalize_gene_symbol("", None, i))
            .collect::<Vec<_>>()
    };

    if gene_ids.len() != n_genes {
        if strict {
            return Err(ScioError::new(
                ErrorCode::DimensionMismatch,
                format!(
                    "gene vector length {} != matrix n_genes {}",
                    gene_ids.len(),
                    n_genes
                ),
            ));
        }
        gene_ids.resize_with(n_genes, || "".to_string());
        for (i, g) in gene_ids.iter_mut().enumerate() {
            if g.is_empty() {
                *g = normalize_gene_symbol("", None, i);
            }
        }
    }

    let mut barcodes = if let Some(path) = ds.barcodes.as_ref() {
        parse_barcodes(path)?
    } else {
        (0..n_cells)
            .map(|i| normalize_barcode("", i))
            .collect::<Vec<_>>()
    };

    if barcodes.len() != n_cells {
        if strict {
            return Err(ScioError::new(
                ErrorCode::DimensionMismatch,
                format!(
                    "barcode vector length {} != matrix n_cells {}",
                    barcodes.len(),
                    n_cells
                ),
            ));
        }
        barcodes.resize_with(n_cells, || "".to_string());
        for (i, b) in barcodes.iter_mut().enumerate() {
            if b.is_empty() {
                *b = normalize_barcode("", i);
            }
        }
    }

    let mut col_ptr = Vec::with_capacity(n_cells + 1);
    let mut row_idx = Vec::with_capacity(stats.nnz);
    let mut values = Vec::with_capacity(stats.nnz);
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
        n_cells,
        n_genes,
        col_ptr,
        row_idx,
        values,
    };
    matrix.validate()?;

    let metadata = InputMetadata {
        format: "mtx10x".to_string(),
        n_cells,
        n_genes,
        gene_ids: gene_ids.clone(),
        gene_symbols: gene_ids,
        barcodes,
        stats,
    };

    Ok((metadata, matrix))
}

fn parse_matrix_market(
    reader: &mut BufReader<Box<dyn Read>>,
    source: &Path,
    strict: bool,
) -> ScioResult<(usize, usize, Vec<Vec<(usize, f32)>>, MatrixStats)> {
    let mut n_rows = None::<usize>;
    let mut n_cols = None::<usize>;
    let mut columns: Vec<Vec<(usize, f32)>> = Vec::new();

    let mut nnz = 0usize;
    let mut total = 0f64;
    let mut min = f32::MAX;
    let mut max = f32::MIN;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }

        if n_rows.is_none() {
            let header = t.split_whitespace().collect::<Vec<_>>();
            if header.len() != 3 {
                return Err(ScioError::new(
                    ErrorCode::ParseError,
                    format!(
                        "{}: invalid MTX header line {}",
                        source.display(),
                        line_no + 1
                    ),
                ));
            }
            let r = header[0]
                .parse::<usize>()
                .map_err(|_| ScioError::new(ErrorCode::ParseError, "invalid n_rows"))?;
            let c = header[1]
                .parse::<usize>()
                .map_err(|_| ScioError::new(ErrorCode::ParseError, "invalid n_cols"))?;
            n_rows = Some(r);
            n_cols = Some(c);
            columns = vec![Vec::new(); c];
            continue;
        }

        let parts = t.split_whitespace().collect::<Vec<_>>();
        if parts.len() < 3 {
            return Err(ScioError::new(
                ErrorCode::ParseError,
                format!(
                    "{}: malformed coordinate at line {}",
                    source.display(),
                    line_no + 1
                ),
            ));
        }
        let row_1 = parts[0]
            .parse::<usize>()
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "invalid row index"))?;
        let col_1 = parts[1]
            .parse::<usize>()
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "invalid col index"))?;
        let val = parts[2]
            .parse::<f32>()
            .map_err(|_| ScioError::new(ErrorCode::ParseError, "invalid value"))?;

        if row_1 == 0 || col_1 == 0 {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                "MTX is 1-based; found zero index",
            ));
        }

        let row = row_1 - 1;
        let col = col_1 - 1;
        if row >= n_rows.unwrap_or(0) || col >= n_cols.unwrap_or(0) {
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!("index out of range at line {}", line_no + 1),
                ));
            }
            continue;
        }

        if val != 0.0 {
            columns[col].push((row, val));
            nnz += 1;
            total += val as f64;
            min = min.min(val);
            max = max.max(val);
        }
    }

    let rows = n_rows.ok_or_else(|| ScioError::new(ErrorCode::ParseError, "missing MTX header"))?;
    let cols = n_cols.ok_or_else(|| ScioError::new(ErrorCode::ParseError, "missing MTX header"))?;

    let sparsity = 1.0 - (nnz as f64 / ((rows * cols) as f64));
    let stats = MatrixStats {
        nnz,
        total_counts: total,
        min_count: if nnz > 0 { min } else { 0.0 },
        max_count: if nnz > 0 { max } else { 0.0 },
        sparsity,
    };

    Ok((rows, cols, columns, stats))
}

fn parse_features(path: &Path) -> ScioResult<Vec<String>> {
    let reader = open_maybe_gz_existing(path)?;
    let mut out = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let cols = t.split('\t').collect::<Vec<_>>();
        let id = cols.first().copied().unwrap_or_default();
        let symbol = cols.get(1).copied();
        out.push(normalize_gene_symbol(id, symbol, line_no));
    }
    Ok(out)
}

fn parse_legacy_genes(path: &Path) -> ScioResult<Vec<String>> {
    let reader = open_maybe_gz_existing(path)?;
    let mut out = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let cols = t.split('\t').collect::<Vec<_>>();
        let id = cols.first().copied().unwrap_or_default();
        let symbol = cols.get(1).copied();
        out.push(normalize_gene_symbol(id, symbol, line_no));
    }
    Ok(out)
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

    let prefix = detect_prefix(&input_dir)?;
    let matrix = candidate_path(&input_dir, prefix.as_deref(), "matrix.mtx");
    let features = candidate_path(&input_dir, prefix.as_deref(), "features.tsv");
    let genes = candidate_path(&input_dir, prefix.as_deref(), "genes.tsv");
    let barcodes = candidate_path(&input_dir, prefix.as_deref(), "barcodes.tsv");

    let matrix = choose_existing(&matrix).ok_or_else(|| {
        ScioError::new(
            ErrorCode::MissingFile,
            format!("missing matrix.mtx(.gz) in {}", input_dir.display()),
        )
    })?;

    Ok(MtxDatasetPaths {
        input_dir,
        prefix,
        matrix,
        features: choose_existing(&features),
        genes: choose_existing(&genes),
        barcodes: choose_existing(&barcodes),
    })
}

pub fn detect_prefix(input_dir: &Path) -> ScioResult<Option<String>> {
    let mut prefixes = std::collections::BTreeSet::new();
    for entry in std::fs::read_dir(input_dir)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else { continue };
        if let Some(p) = extract_prefix(name) {
            prefixes.insert(p.to_string());
        }
    }
    if prefixes.len() > 1 {
        return Err(ScioError::new(
            ErrorCode::ValidationError,
            format!("multiple dataset prefixes in {}", input_dir.display()),
        ));
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
    let suffixes = [
        "_matrix.mtx",
        "_matrix.mtx.gz",
        ".matrix.mtx",
        ".matrix.mtx.gz",
        "_features.tsv",
        "_features.tsv.gz",
        ".features.tsv",
        ".features.tsv.gz",
        "_barcodes.tsv",
        "_barcodes.tsv.gz",
        ".barcodes.tsv",
        ".barcodes.tsv.gz",
        "_genes.tsv",
        "_genes.tsv.gz",
        ".genes.tsv",
        ".genes.tsv.gz",
    ];
    for suffix in suffixes {
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
        Some(p) => {
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
        None => input_dir.join(name),
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
    })?;
    if existing
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        let f = File::open(&existing)?;
        Ok(BufReader::new(Box::new(GzDecoder::new(f))))
    } else {
        let f = File::open(&existing)?;
        Ok(BufReader::new(Box::new(f)))
    }
}
