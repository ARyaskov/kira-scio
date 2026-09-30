//! 10x Matrix Market reader (`matrix.mtx[.gz]` + features/genes + barcodes).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use rustc_hash::FxHashSet;
use tracing::warn;

use crate::error::{ErrorCode, ScioError, ScioResult};
use crate::model::{
    CountMismatch, IngestReport, InputMetadata, MatrixKind, MatrixStats, Provenance, ShapeProbe,
    SoaCscMatrix,
};
use crate::normalize::{normalize_barcode, normalize_gene_id, normalize_gene_symbol, strip_bom};

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
    let mut parsed = parse_matrix_market(matrix_reader, &ds.matrix, strict)?;
    let mut report = std::mem::take(&mut parsed.report);

    let (mut gene_ids, mut gene_symbols, feature_types) =
        if let Some(features_path) = ds.features.as_ref() {
            parse_features(features_path, strict, &mut report)?
        } else if let Some(genes_path) = ds.genes.as_ref() {
            parse_features(genes_path, strict, &mut report)?
        } else {
            let synth: Vec<String> = (0..parsed.n_genes)
                .map(|i| normalize_gene_id("", None, i))
                .collect();
            (synth.clone(), synth, None)
        };
    let mut barcodes = if let Some(path) = ds.barcodes.as_ref() {
        parse_barcodes(path, &mut report)?
    } else {
        (0..parsed.n_cells).map(normalize_barcode_idx).collect()
    };

    // Some exporters write the matrix as cells x genes. When both label
    // files are present and their lengths match the swapped dimensions but
    // not the declared ones, the orientation is unambiguous: transpose.
    let both_label_files = (ds.features.is_some() || ds.genes.is_some()) && ds.barcodes.is_some();
    let declared_fit = gene_ids.len() == parsed.n_genes && barcodes.len() == parsed.n_cells;
    let swapped_fit = gene_ids.len() == parsed.n_cells
        && barcodes.len() == parsed.n_genes
        && parsed.n_genes != parsed.n_cells;
    if both_label_files && !declared_fit && swapped_fit {
        parsed.transpose();
        report.transposed = true;
    }

    // A type column shorter than the matrix cannot be trusted for filtering.
    let feature_types = feature_types.filter(|t| t.len() == parsed.n_genes);

    report.relabeled_genes = fix_length(
        &mut gene_ids,
        parsed.n_genes,
        strict,
        "gene",
        &ds.matrix,
        |i| normalize_gene_id("", None, i),
    )?;
    // Symbols come from the same file as ids, so they share the mismatch.
    fix_length(
        &mut gene_symbols,
        parsed.n_genes,
        strict,
        "gene_symbol",
        &ds.matrix,
        |i| normalize_gene_symbol("", None, i),
    )?;
    report.relabeled_barcodes = fix_length(
        &mut barcodes,
        parsed.n_cells,
        strict,
        "barcode",
        &ds.matrix,
        normalize_barcode_idx,
    )?;

    let (matrix, merged_duplicates) = parsed.into_csc();
    report.merged_duplicates = merged_duplicates;
    matrix.validate()?;
    let stats = MatrixStats::from_matrix(&matrix);
    log_report(&ds.matrix, &report);

    let metadata = InputMetadata {
        format: "mtx10x".to_string(),
        n_cells: matrix.n_cells,
        n_genes: matrix.n_genes,
        gene_ids,
        gene_symbols,
        barcodes,
        stats,
        feature_types,
        marginals: Default::default(),
        report,
        provenance: Provenance {
            dialect: if ds.features.is_some() {
                "mtx-v3"
            } else if ds.genes.is_some() {
                "mtx-v2"
            } else {
                "mtx"
            }
            .to_string(),
            dataset_prefix: ds.prefix.clone(),
            matrix_kind: MatrixKind::from_path(&ds.matrix),
            source_path: ds.matrix,
            matrix_source: None,
        },
    };

    Ok((metadata, matrix))
}

/// Dimensions from the Matrix Market header plus the label files, applying
/// the same transposition rule as [`read_mtx`], without reading entries.
pub(crate) fn read_shape(path: &Path) -> ScioResult<ShapeProbe> {
    let ds = discover(path)?;
    let reader = open_maybe_gz_existing(&ds.matrix)?;
    let mut n_rows = None::<usize>;
    let mut n_cols = None::<usize>;
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = if line_no == 0 {
            strip_bom(&line)
        } else {
            &line
        };
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }
        let mut header = t.split_whitespace();
        let parse = |tok: Option<&str>, what: &str| -> ScioResult<usize> {
            tok.and_then(|s| s.parse::<usize>().ok()).ok_or_else(|| {
                ScioError::new(ErrorCode::ParseError, format!("invalid {what}"))
                    .with_path(ds.matrix.clone())
            })
        };
        n_rows = Some(parse(header.next(), "n_rows")?);
        n_cols = Some(parse(header.next(), "n_cols")?);
        break;
    }
    let (Some(mut n_genes), Some(mut n_cells)) = (n_rows, n_cols) else {
        return Err(ScioError::new(ErrorCode::ParseError, "missing MTX header")
            .with_path(ds.matrix.clone()));
    };

    let mut scratch = IngestReport::default();
    let features_path = ds.features.as_ref().or(ds.genes.as_ref());
    let (n_feature_labels, feature_types) = match features_path {
        Some(p) => {
            let (ids, _, types) = parse_features(p, false, &mut scratch)?;
            (Some(ids.len()), types)
        }
        None => (None, None),
    };
    let n_barcode_labels = match ds.barcodes.as_ref() {
        Some(p) => Some(parse_barcodes(p, &mut scratch)?.len()),
        None => None,
    };
    if let (Some(nf), Some(nb)) = (n_feature_labels, n_barcode_labels) {
        let declared_fit = nf == n_genes && nb == n_cells;
        let swapped_fit = nf == n_cells && nb == n_genes && n_genes != n_cells;
        if !declared_fit && swapped_fit {
            std::mem::swap(&mut n_genes, &mut n_cells);
        }
    }
    let feature_types = feature_types.filter(|t| t.len() == n_genes);
    Ok(ShapeProbe {
        n_cells,
        n_genes,
        feature_types,
    })
}

fn normalize_barcode_idx(idx: usize) -> String {
    crate::normalize::synth_barcode(idx)
}

/// One warning per repair category, with counts rather than per-entry noise.
pub(crate) fn log_report(source: &Path, report: &IngestReport) {
    let path = source.display();
    if report.dropped_out_of_range > 0 {
        warn!(%path, count = report.dropped_out_of_range, "dropped entries with out-of-range coordinates");
    }
    if report.dropped_non_finite > 0 {
        warn!(%path, count = report.dropped_non_finite, "dropped non-finite values");
    }
    if report.merged_duplicates > 0 {
        warn!(%path, count = report.merged_duplicates, "duplicate coordinates were summed");
    }
    if let Some(m) = report.entry_count_mismatch {
        warn!(%path, expected = m.expected, found = m.found, "entry count does not match header");
    }
    if let Some(m) = report.relabeled_genes {
        warn!(%path, expected = m.expected, found = m.found, "gene label count resized to matrix");
    }
    if let Some(m) = report.relabeled_barcodes {
        warn!(%path, expected = m.expected, found = m.found, "barcode count resized to matrix");
    }
    if report.transposed {
        warn!(%path, "matrix was stored cells x genes; transposed to genes x cells");
    }
}

/// Resizes a label vector to the matrix dimension. Returns the recorded
/// mismatch in lenient mode; strict mode errors instead.
pub(crate) fn fix_length(
    out: &mut Vec<String>,
    expected: usize,
    strict: bool,
    label: &'static str,
    source: &Path,
    synth: impl Fn(usize) -> String,
) -> ScioResult<Option<CountMismatch>> {
    if out.len() == expected {
        return Ok(None);
    }
    let mismatch = CountMismatch {
        expected,
        found: out.len(),
    };
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
    out.resize_with(expected, String::new);
    for (i, v) in out.iter_mut().enumerate() {
        if v.is_empty() {
            *v = synth(i);
        }
    }
    Ok(Some(mismatch))
}

/// Intermediate triplet form; `into_csc()` canonicalizes into CSC.
struct ParsedMtx {
    n_genes: usize,
    n_cells: usize,
    /// `(col, row, value)` triplets in file order.
    triplets: Vec<(u32, u32, f32)>,
    /// Repairs recorded while scanning the matrix file.
    report: IngestReport,
}

impl ParsedMtx {
    /// Swaps the gene and cell axes in place.
    fn transpose(&mut self) {
        std::mem::swap(&mut self.n_genes, &mut self.n_cells);
        for t in &mut self.triplets {
            std::mem::swap(&mut t.0, &mut t.1);
        }
    }

    /// Returns the canonical matrix and the number of merged duplicate
    /// coordinates (summed per Matrix Market convention).
    fn into_csc(self) -> (SoaCscMatrix, usize) {
        SoaCscMatrix::from_triplets(self.n_cells, self.n_genes, self.triplets)
    }
}

/// Value field declared by the Matrix Market banner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MtxField {
    Integer,
    Real,
    /// Presence-only entries: two tokens per line, value taken as 1.
    Pattern,
}

/// Parses `%%MatrixMarket matrix coordinate <field> <symmetry>`.
fn parse_banner(line: &str, source: &Path) -> ScioResult<MtxField> {
    let unsupported = |msg: String| {
        ScioError::new(ErrorCode::UnsupportedFormat, msg).with_path(source.to_path_buf())
    };
    let mut it = line.split_whitespace().skip(1).map(str::to_ascii_lowercase);
    let object = it.next().unwrap_or_default();
    let format = it.next().unwrap_or_default();
    let field = it.next().unwrap_or_default();
    let symmetry = it.next().unwrap_or_else(|| "general".to_string());
    if object != "matrix" {
        return Err(unsupported(format!(
            "Matrix Market object `{object}` is not supported (expected `matrix`)"
        )));
    }
    if format != "coordinate" {
        return Err(unsupported(format!(
            "Matrix Market format `{format}` is not supported (expected `coordinate`; \
             dense `array` files must be converted first)"
        )));
    }
    if symmetry != "general" {
        return Err(unsupported(format!(
            "Matrix Market symmetry `{symmetry}` is not supported (expected `general`)"
        )));
    }
    match field.as_str() {
        "integer" => Ok(MtxField::Integer),
        "real" | "double" => Ok(MtxField::Real),
        "pattern" => Ok(MtxField::Pattern),
        other => Err(unsupported(format!(
            "Matrix Market field `{other}` is not supported (expected integer, real or pattern)"
        ))),
    }
}

/// Advances past ASCII whitespace and returns the next token, if any.
fn take_token<'a>(cur: &mut &'a [u8]) -> Option<&'a [u8]> {
    let start = cur.iter().position(|b| !b.is_ascii_whitespace())?;
    let rest = &cur[start..];
    let end = rest
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(rest.len());
    *cur = &rest[end..];
    Some(&rest[..end])
}

/// Decimal unsigned integer with overflow checking.
fn parse_uint(tok: &[u8]) -> Option<usize> {
    if tok.is_empty() {
        return None;
    }
    let mut v: usize = 0;
    for &b in tok {
        if !b.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((b - b'0') as usize)?;
    }
    Some(v)
}

/// Optionally signed decimal integer as `f32`; `None` when the token is not
/// a plain integer (caller falls back to float parsing).
fn parse_int_value(tok: &[u8]) -> Option<f32> {
    let (neg, digits) = match tok.first()? {
        b'-' => (true, &tok[1..]),
        b'+' => (false, &tok[1..]),
        _ => (false, tok),
    };
    let magnitude = parse_uint(digits)? as f32;
    Some(if neg { -magnitude } else { magnitude })
}

fn parse_value(tok: &[u8], field: MtxField) -> Option<f32> {
    if field == MtxField::Integer
        && let Some(v) = parse_int_value(tok)
    {
        return Some(v);
    }
    std::str::from_utf8(tok).ok()?.parse::<f32>().ok()
}

fn parse_matrix_market(
    mut reader: BufReader<Box<dyn Read>>,
    source: &Path,
    strict: bool,
) -> ScioResult<ParsedMtx> {
    let mut n_rows = None::<usize>;
    let mut n_cols = None::<usize>;
    let mut nnz_hint = None::<usize>;
    let mut entries_seen = 0usize;
    let mut triplets: Vec<(u32, u32, f32)> = Vec::new();
    let mut report = IngestReport::default();
    let mut field: Option<MtxField> = None;

    let malformed = |what: &str, line_no: usize| {
        ScioError::new(
            ErrorCode::ParseError,
            format!("malformed {what} at line {line_no}"),
        )
        .with_path(source.to_path_buf())
    };

    // One reusable buffer instead of a String per line.
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    let mut line_no = 0usize;
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        line_no += 1;
        let mut line: &[u8] = &buf;
        if line_no == 1 && line.starts_with(UTF8_BOM_BYTES) {
            report.bom_stripped = true;
            line = &line[UTF8_BOM_BYTES.len()..];
        }
        let line = line.trim_ascii();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'%' {
            if field.is_none() && n_rows.is_none() {
                let text = String::from_utf8_lossy(line);
                if text
                    .get(..14)
                    .is_some_and(|p| p.eq_ignore_ascii_case("%%MatrixMarket"))
                {
                    field = Some(parse_banner(&text, source)?);
                }
            }
            continue;
        }

        if n_rows.is_none() {
            if field.is_none() {
                if strict {
                    return Err(ScioError::new(
                        ErrorCode::ParseError,
                        "missing %%MatrixMarket banner (use strict=false to assume `real general`)",
                    )
                    .with_path(source.to_path_buf()));
                }
                field = Some(MtxField::Real);
            }
            let mut cur = line;
            let r = take_token(&mut cur)
                .and_then(parse_uint)
                .ok_or_else(|| header_err(source, line_no - 1))?;
            let c = take_token(&mut cur)
                .and_then(parse_uint)
                .ok_or_else(|| header_err(source, line_no - 1))?;
            // Third header field: number of entries that follow. Used to
            // detect truncated files and to pre-size the triplet buffer.
            match take_token(&mut cur) {
                Some(tok) => {
                    let hint = parse_uint(tok).ok_or_else(|| {
                        ScioError::new(ErrorCode::ParseError, "invalid nnz in MTX header")
                            .with_path(source.to_path_buf())
                    })?;
                    nnz_hint = Some(hint);
                    triplets.reserve(hint.min(MAX_RESERVED_ENTRIES));
                }
                None if strict => return Err(header_err(source, line_no - 1)),
                None => {}
            }
            SoaCscMatrix::check_dims(c, r).map_err(|e| e.with_path(source.to_path_buf()))?;
            n_rows = Some(r);
            n_cols = Some(c);
            continue;
        }
        let field = field.unwrap_or(MtxField::Real);

        entries_seen += 1;
        let mut cur = line;
        let row_1 = take_token(&mut cur)
            .and_then(parse_uint)
            .ok_or_else(|| malformed("coordinate", line_no))?;
        let col_1 = take_token(&mut cur)
            .and_then(parse_uint)
            .ok_or_else(|| malformed("coordinate", line_no))?;
        let val = match field {
            MtxField::Pattern => 1.0,
            f => take_token(&mut cur)
                .and_then(|tok| parse_value(tok, f))
                .ok_or_else(|| malformed("value", line_no))?,
        };

        if row_1 == 0 || col_1 == 0 {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!("MTX is 1-based; found zero index at line {line_no}"),
            )
            .with_path(source.to_path_buf()));
        }

        let row = row_1 - 1;
        let col = col_1 - 1;
        if row >= n_rows.unwrap_or(0) || col >= n_cols.unwrap_or(0) {
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!("index out of range at line {line_no} (use strict=false to drop)"),
                )
                .with_path(source.to_path_buf()));
            }
            report.dropped_out_of_range += 1;
            continue;
        }

        if !val.is_finite() {
            if strict {
                return Err(ScioError::new(
                    ErrorCode::ValidationError,
                    format!("non-finite value at line {line_no} (use strict=false to ignore)"),
                )
                .with_path(source.to_path_buf()));
            }
            report.dropped_non_finite += 1;
            continue;
        }

        if val == 0.0 {
            report.explicit_zeros += 1;
        } else {
            triplets.push((col as u32, row as u32, val));
        }
    }

    if let Some(expected) = nnz_hint
        && expected != entries_seen
    {
        if strict {
            return Err(ScioError::new(
                ErrorCode::ParseError,
                format!(
                    "MTX header declares {expected} entries but {entries_seen} were found \
                     (truncated or corrupt file; use strict=false to accept)"
                ),
            )
            .with_path(source.to_path_buf()));
        }
        report.entry_count_mismatch = Some(CountMismatch {
            expected,
            found: entries_seen,
        });
    }

    let rows = n_rows.ok_or_else(|| {
        ScioError::new(ErrorCode::ParseError, "missing MTX header").with_path(source.to_path_buf())
    })?;
    let cols = n_cols.ok_or_else(|| {
        ScioError::new(ErrorCode::ParseError, "missing MTX header").with_path(source.to_path_buf())
    })?;

    Ok(ParsedMtx {
        n_genes: rows,
        n_cells: cols,
        triplets,
        report,
    })
}

const UTF8_BOM_BYTES: &[u8] = "\u{feff}".as_bytes();

/// Upper bound on entries pre-allocated from the header hint so a bogus
/// header cannot trigger a multi-gigabyte allocation up front.
const MAX_RESERVED_ENTRIES: usize = 1 << 26;

fn header_err(source: &Path, line_no: usize) -> ScioError {
    ScioError::new(
        ErrorCode::ParseError,
        format!("invalid MTX header line {}", line_no + 1),
    )
    .with_path(source.to_path_buf())
}

/// Returns `(gene_ids, gene_symbols, feature_types)`; single-column rows
/// reuse the id. `feature_types` is `Some` when at least one row carries the
/// third (`feature_type`) column of 10x v3 `features.tsv`; rows without it
/// get an empty string.
#[allow(clippy::type_complexity)]
fn parse_features(
    path: &Path,
    strict: bool,
    report: &mut IngestReport,
) -> ScioResult<(Vec<String>, Vec<String>, Option<Vec<String>>)> {
    let reader = open_maybe_gz_existing(path)?;
    let mut ids = Vec::<String>::new();
    let mut symbols = Vec::<String>::new();
    let mut types = Vec::<String>::new();
    let mut saw_type_column = false;
    let mut single_column_row: Option<usize> = None;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = if line_no == 0 {
            let stripped = strip_bom(&line);
            report.bom_stripped |= stripped.len() != line.len();
            stripped
        } else {
            &line
        };
        let t = line.trim_end_matches(['\r', '\n']);
        if t.trim().is_empty() {
            continue;
        }
        let mut it = t.splitn(4, '\t');
        let id = it.next().unwrap_or("").trim();
        let sym = it.next().map(|s| s.trim());
        let ty = it.next().map(|s| s.trim()).unwrap_or("");
        if sym.is_none() {
            single_column_row.get_or_insert(line_no);
        }
        saw_type_column |= !ty.is_empty();
        ids.push(normalize_gene_id(id, sym, ids.len()));
        symbols.push(normalize_gene_symbol(id, sym, symbols.len()));
        types.push(ty.to_string());
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

    Ok((ids, symbols, saw_type_column.then_some(types)))
}

fn parse_barcodes(path: &Path, report: &mut IngestReport) -> ScioResult<Vec<String>> {
    let reader = open_maybe_gz_existing(path)?;
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = if i == 0 {
            let stripped = strip_bom(&line);
            report.bom_stripped |= stripped.len() != line.len();
            stripped
        } else {
            &line
        };
        if line.trim().is_empty() {
            continue;
        }
        out.push(normalize_barcode(line, i));
    }
    Ok(out)
}

/// Locates the MTX triplet for `path`.
///
/// `path` may be the dataset directory or the `matrix.mtx[.gz]` file itself.
/// Given a file, its own name fixes the dataset prefix, so a directory that
/// holds several prefixed datasets is not ambiguous. Given a directory, the
/// files must agree on a single prefix.
pub fn discover(path: &Path) -> ScioResult<MtxDatasetPaths> {
    let (input_dir, explicit_prefix) = if path.is_dir() {
        (path.to_path_buf(), None)
    } else {
        let dir = path
            .parent()
            .ok_or_else(|| ScioError::new(ErrorCode::InvalidInputPath, "input has no parent"))?
            .to_path_buf();
        let file_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        // `Some(Some(p))` prefixed, `Some(None)` unprefixed file.
        (dir, Some(extract_prefix(file_name).map(str::to_string)))
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
    let prefix = match explicit_prefix {
        Some(p) => p,
        None if prefixes.len() > 1 => {
            return Err(ScioError::new(
                ErrorCode::ValidationError,
                format!(
                    "multiple dataset prefixes in {} ({}); pass the matrix file explicitly",
                    input_dir.display(),
                    prefixes.iter().cloned().collect::<Vec<_>>().join(", ")
                ),
            )
            .with_path(input_dir));
        }
        None => prefixes.into_iter().next(),
    };

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
