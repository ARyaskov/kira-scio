# kira-scio

`kira-scio` is a standalone **library crate** for deterministic single-cell input ingestion in the Kira stack.

## Scope

- Auto-detect input format:
  - 10x Genomics MTX (v2/v3)
  - BD Rhapsody WTA (`raw_counts.tsv(.gz)`, `*_raw_counts.tsv(.gz)`)
  - Dense TSV/CSV (`.tsv/.csv`, gz supported)
  - H5AD (feature-gated)
  - loom (optional placeholder)
- Unified canonical model:
  - cell × gene
  - sparse CSC (SoA layout)
  - deterministic ordering
- Unified API:
  - `read_metadata()`
  - `read_matrix()`
  - `read_all()`
- Strict error taxonomy with stable error codes.
- Zero biology in this crate.

## API

```rust
use kira_scio::Reader;

let reader = Reader::new("/path/to/input");
let metadata = reader.read_metadata()?;
let matrix = reader.read_matrix()?;
let all = reader.read_all()?;
```

## Notes

- Parsing is streaming-oriented for text formats (line-by-line).
- `.gz` variants are supported where applicable.
- MTX prefix variants are supported for both underscore and dot naming styles.
- H5AD support requires `--features h5ad`.

