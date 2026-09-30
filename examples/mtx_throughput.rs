//! Reproducible throughput check for the MTX reader.
//!
//! Generates a synthetic 10x-style dataset (default 5 000 cells x 30 000
//! genes, 2 000 entries per cell, written column-major like Cell Ranger),
//! then times `read_shape`, `read_all` and a second `read_metadata` served
//! from the cache. Numbers are wall-clock and noisy; compare medians of a
//! few runs on the same machine.
//!
//! ```bash
//! cargo run --release --example mtx_throughput -- [cells] [genes] [entries_per_cell]
//! ```

use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use kira_scio::Reader;

fn generate(dir: &Path, cells: usize, genes: usize, per_cell: usize) -> u64 {
    let mut w = BufWriter::new(std::fs::File::create(dir.join("matrix.mtx")).unwrap());
    writeln!(
        w,
        "%%MatrixMarket matrix coordinate integer general\n{genes} {cells} {}",
        cells * per_cell
    )
    .unwrap();
    // xorshift keeps the data deterministic across runs.
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let stride = genes / per_cell;
    for c in 1..=cells {
        for k in 0..per_cell {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let g = k * stride + (x as usize % stride) + 1;
            writeln!(w, "{g} {c} {}", 1 + x % 20).unwrap();
        }
    }
    w.flush().unwrap();
    let mut f = BufWriter::new(std::fs::File::create(dir.join("features.tsv")).unwrap());
    for g in 0..genes {
        writeln!(f, "ENSG{g:011}\tGENE{g}\tGene Expression").unwrap();
    }
    let mut b = BufWriter::new(std::fs::File::create(dir.join("barcodes.tsv")).unwrap());
    for c in 0..cells {
        writeln!(b, "BC{c:016}-1").unwrap();
    }
    std::fs::metadata(dir.join("matrix.mtx")).unwrap().len()
}

fn main() {
    let args: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("numeric argument"))
        .collect();
    let cells = args.first().copied().unwrap_or(5_000);
    let genes = args.get(1).copied().unwrap_or(30_000);
    let per_cell = args.get(2).copied().unwrap_or(2_000).min(genes);

    let dir = std::env::temp_dir().join(format!("kira_scio_throughput_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();

    let t = Instant::now();
    let bytes = generate(&dir, cells, genes, per_cell);
    println!(
        "generated {cells} x {genes}, {} entries, {:.0} MB in {:.2}s",
        cells * per_cell,
        bytes as f64 / 1e6,
        t.elapsed().as_secs_f64()
    );

    let reader = Reader::new(&dir);
    let t = Instant::now();
    let shape = reader.read_shape().unwrap();
    println!(
        "read_shape:    {:?} in {:.3}s",
        shape,
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let data = reader.read_all().unwrap();
    let el = t.elapsed().as_secs_f64();
    println!(
        "read_all:      {:.2}s  ({:.0} MB/s, {:.1} M entries/s), nnz={}, sparsity={:.4}",
        el,
        bytes as f64 / 1e6 / el,
        data.metadata.stats.nnz as f64 / 1e6 / el,
        data.metadata.stats.nnz,
        data.metadata.stats.sparsity
    );

    let t = Instant::now();
    let md = reader.read_metadata().unwrap();
    println!(
        "read_metadata: {:.4}s from cache ({} cells)",
        t.elapsed().as_secs_f64(),
        md.n_cells
    );

    std::fs::remove_dir_all(&dir).ok();
}
