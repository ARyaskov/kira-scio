# kira-scio

`kira-scio` is a standalone **library crate** for deterministic single-cell input ingestion in the Kira stack.

## Scope

- Auto-detect input format:
  - 10x Genomics MTX (v2/v3), including prefixed (`<sample>_matrix.mtx.gz`, `<sample>.matrix.mtx`) triplets
  - BD Rhapsody WTA: Sequence Analysis Pipeline `*_{DBEC,RSEC}_MolsPerCell.csv(.gz)` (cell-major, `Cell_Index` first column) and the legacy `raw_counts.tsv(.gz)` convention
  - Dense TSV/CSV (`.tsv/.csv`, gz supported), gene-major or cell-major by header label
  - H5AD (feature `h5ad`)
  - 10x Cell Ranger HDF5 `*_feature_bc_matrix.h5`, v2 and v3 layouts (feature `tenx-h5`)
  - loom (optional placeholder)
- Unified canonical model:
  - cell × gene
  - sparse CSC (SoA layout), canonical form: rows strictly increasing within a column, duplicate coordinates summed, no explicit zeros
  - deterministic ordering
- Unified API:
  - `read_shape()` — `(n_cells, n_genes)` from headers and label files, no entries parsed
  - `read_metadata()` — labels, statistics, marginals, report, provenance (cached after a full parse)
  - `read_matrix()`
  - `read_all()`
- Strict error taxonomy with stable error codes.
- Zero biology in this crate.

## API

```rust
use kira_scio::Reader;

let reader = Reader::new("/path/to/input");
let (n_cells, n_genes) = reader.read_shape()?;
let metadata = reader.read_metadata()?;
let matrix = reader.read_matrix()?;
let all = reader.read_all()?;

// What the reader had to tolerate or repair.
let report = &all.metadata.report;
assert!(report.is_lossless());

// Whole-matrix summary and per-cell / per-gene marginals.
assert!(all.metadata.stats.is_integer, "not raw counts");
let library_size = &all.metadata.marginals.cell_total_counts;
let genes_detected = &all.metadata.marginals.cell_n_features;
```

## Options

```rust
use kira_scio::{FeatureTypeFilter, H5adSource, Reader, ReaderOptions};

let reader = Reader::with_options("/path/to/input", ReaderOptions {
    strict: true,
    force_format: None,
    // Keep only `Gene Expression` rows of a Feature Barcoding matrix.
    feature_types: FeatureTypeFilter::GeneExpression,
    // For h5ad: read raw counts from `raw/X` instead of a normalized `X`.
    h5ad_source: H5adSource::RawX,
    // Keep Ensembl ids exactly as written (default strips `.N` versions).
    strip_ensembl_versions: false,
    // scanpy-style `-1`, `-2` suffixes for repeated gene labels.
    make_gene_labels_unique: true,
});
```

- `feature_types`: `All` (default), `GeneExpression`, or `Only([...])`. Applies when the source declares feature types (10x v3 `features.tsv`, 10x `.h5`, anndata `var/feature_types`); the excluded count is recorded in `report.excluded_features`.
- `h5ad_source`: `X` (default), `RawX` (labels from `raw/var`), or `Layer(name)`. A missing source is a `MissingFile` error.
- `strip_ensembl_versions` (default `true`): strips the `.N` suffix from Ensembl stable ids of any species and feature type (`ENSG`, `ENSMUSG`, `ENSDARG`, `ENST`, …).
- `make_gene_labels_unique` (default `false`): renames repeated gene ids and symbols `X`, `X-1`, `X-2`, … The original duplicates stay listed in `report.duplicate_gene_ids`; duplicates are detected the same way for every format.

## Provenance

`InputMetadata::provenance` records the resolved source file, a dialect tag (`mtx-v3`, `mtx-v2`, `bd-molspercell-dbec`, `dense-cell-major`, `anndata-0.1.0`, `cellranger-h5-v3`, …), the MEX dataset prefix, the AnnData matrix source, and `MatrixKind::{Raw, Filtered, Unknown}` inferred from Cell Ranger / BD naming (`raw_feature_bc_matrix`, `filtered_…`, `_Unfiltered_`). Empty-droplet and ambient-RNA methods need `Raw`.

## Statistics

`InputMetadata::stats` (`MatrixStats`) summarizes stored entries: `nnz`, `total_counts`, `min_count`, `max_count`, `sparsity`, plus two integrity flags, `is_integer` and `has_negative`, that tell whether the matrix can be treated as raw counts.

`InputMetadata::marginals` (`Marginals`) holds the first-pass QC covariates computed in one pass over the canonical matrix: per cell `cell_total_counts`, `cell_n_features`, `cell_top_fraction` (share of the 50 largest entries); per gene `gene_total_counts`, `gene_n_cells`. `Marginals::fraction_in_gene_set` computes the per-cell share of any caller-defined gene set (mitochondrial, ribosomal, spike-ins) so the crate stays free of biology.

## Strict and lenient mode

`Reader::new` is strict. `Reader::with_options` with `strict: false` is lenient. The rules are the same for every format and are recorded in `InputMetadata::report` (`IngestReport`):

| Situation | Strict | Lenient |
|---|---|---|
| Coordinate outside the declared shape | error | dropped, counted |
| NaN / infinite value | error | dropped, counted |
| MTX header entry count ≠ entries found (truncation) | error | accepted, recorded |
| Label vector length ≠ matrix dimension | error | resized with synthesized labels, recorded |
| Duplicate coordinates | summed, counted | summed, counted |
| Explicit zeros in a sparse input | dropped, counted | dropped, counted |
| Duplicate gene ids / barcodes | kept as separate rows/columns, listed | same |
| UTF-8 BOM on a text file | stripped, flagged | same |
| MTX stored cells x genes, both label files fit the swapped shape | transposed, flagged | same |
| Structural corruption (inconsistent `indptr`, unreadable header) | error | error |

## HDF5-based formats

Build with `--features h5ad` and/or `--features tenx-h5`. The bindings are `hdf5-metno` (HDF5 1.8 through 2.x); point `HDF5_DIR` at your installation if it is not on the default search path (for example `HDF5_DIR=/opt/homebrew/opt/hdf5`).

### 10x `.h5`

Cell Ranger ≥ 3 `/matrix` layout (`features/{id,name,feature_type}`) and the Cell Ranger 2 single-genome layout (`genes`, `gene_names`). Multi-genome Cell Ranger 2 files are rejected. Fixtures: `tests/fixtures/tenx_h5` (`generate.py`, requires `h5py`).

### H5AD

Supported layouts:

- `/X` as `csr_matrix`, `csc_matrix` (also the legacy `h5sparse_format` attribute) or a dense 2-D dataset; int32/int64 index arrays; integer or float data.
- `obs`/`var` index resolved through the group's `_index` attribute, with conventional names as fallback.
- String columns as variable- or fixed-length string datasets, `nullable-string-array` groups (anndata ≥ 0.10) or `categorical` groups.
- Gene symbols from the first of `gene_symbols`, `feature_name`, `gene_symbol`, `gene_name`, `symbol`.

`ReaderOptions::h5ad_source` selects `/X`, `raw/X` or `layers/<name>`. `.h5ad.gz` must be decompressed first.

Test fixtures live in `tests/fixtures/h5ad` and are regenerated with `python3 tests/fixtures/h5ad/generate.py tests/fixtures/h5ad` (requires `anndata`, `h5py`, `scipy`).

## Matrix Market

The banner is validated: `matrix coordinate {integer|real|pattern} general` is accepted (`pattern` entries count as 1); `array`, `complex` and non-`general` symmetry are `UnsupportedFormat`. A missing banner is an error in strict mode and means `real general` in lenient mode. Entries are parsed on bytes with a reused buffer; the canonical sort is skipped when entries already arrive column-major, as Cell Ranger and SciPy write them.

A directory must hold a single prefixed dataset; to read one of several (`A_matrix.mtx`, `B_matrix.mtx`), pass the matrix file itself.

## Format detection

Extension first (`.mtx`, `.h5ad`, `.h5`, `.tsv`, `.csv`, `.gz` variants, BD `*_MolsPerCell.csv` and `raw_counts.tsv` names), then content sniffing for other files: a Matrix Market header, a `#`-commented table or a `Cell_Index` first column identify MTX and BD Rhapsody; everything else is a dense table. The value type is not a signal: fractional values only mean the table was normalized upstream, which `MatrixStats::is_integer` reports.

## Development

```bash
cargo test                       # text formats
cargo test --all-features        # + h5ad and 10x .h5 (needs libhdf5; set HDF5_DIR if not on the default path)
cargo clippy --all-features --all-targets -- -D warnings
cargo run --release --example mtx_throughput   # synthetic 10x MTX, prints MB/s for read_all
```

Integration tests write only under a per-test scratch directory that is removed on drop (`tests/common`). Property tests (`tests/proptest_csc.rs`) check the canonical-CSC invariants against a dense reference. CI runs fmt, clippy and the test matrix on Linux and macOS, and verifies that the committed HDF5 fixtures match what `generate.py` produces.

## Notes

- Parsing is streaming-oriented for text formats (line-by-line).
- `.gz` variants are supported where applicable.
- MTX prefix variants are supported for both underscore and dot naming styles.
