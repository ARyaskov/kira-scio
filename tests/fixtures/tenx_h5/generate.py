#!/usr/bin/env python3
"""Regenerates the 10x Genomics HDF5 fixtures used by tests/tenx_h5.rs.

Both files encode the same 3 features x 2 barcodes matrix in the Cell Ranger
layout (CSC, columns = barcodes, rows = features, fixed-length strings):

    feature          AAACCCAAGAAACACT-1  AAACCCAAGAAACCAT-1
    ENSG1 / A                 5                   0
    ENSG2 / B                 0                   3
    CD3_TotalSeqB (ADT)       1                   0

v3_filtered.h5 uses the Cell Ranger >= 3 `matrix/features` group with a
feature_type column; v2_legacy.h5 uses the Cell Ranger 2 per-genome group
with `genes` / `gene_names` datasets (no feature types).

Requires: h5py, numpy.

    python3 tests/fixtures/tenx_h5/generate.py tests/fixtures/tenx_h5
"""
import sys

import h5py
import numpy as np

out = sys.argv[1] if len(sys.argv) > 1 else "."

barcodes = ["AAACCCAAGAAACACT-1", "AAACCCAAGAAACCAT-1"]
ids = ["ENSG1", "ENSG2", "CD3_TotalSeqB"]
names = ["A", "B", "CD3_TotalSeqB"]
feature_types = ["Gene Expression", "Gene Expression", "Antibody Capture"]
# CSC over barcodes: column 0 -> rows 0 and 2; column 1 -> row 1.
data = np.array([5, 1, 3], dtype=np.int32)
indices = np.array([0, 2, 1], dtype=np.int64)
indptr = np.array([0, 2, 3], dtype=np.int64)
shape = np.array([3, 2], dtype=np.int32)


def fixed(strings):
    width = max(len(s) for s in strings)
    return np.array([s.encode() for s in strings], dtype=f"S{width}")


with h5py.File(f"{out}/v3_filtered.h5", "w") as f:
    f.attrs["filetype"] = "matrix"
    f.attrs["version"] = 2
    m = f.create_group("matrix")
    m.create_dataset("barcodes", data=fixed(barcodes))
    m.create_dataset("data", data=data)
    m.create_dataset("indices", data=indices)
    m.create_dataset("indptr", data=indptr)
    m.create_dataset("shape", data=shape)
    feat = m.create_group("features")
    feat.create_dataset("_all_tag_keys", data=fixed(["genome"]))
    feat.create_dataset("id", data=fixed(ids))
    feat.create_dataset("name", data=fixed(names))
    feat.create_dataset("feature_type", data=fixed(feature_types))
    feat.create_dataset("genome", data=fixed(["GRCh38", "GRCh38", ""]))

with h5py.File(f"{out}/v2_legacy.h5", "w") as f:
    g = f.create_group("GRCh38")
    g.create_dataset("barcodes", data=fixed(barcodes))
    g.create_dataset("data", data=data)
    g.create_dataset("indices", data=indices)
    g.create_dataset("indptr", data=indptr)
    g.create_dataset("shape", data=shape)
    g.create_dataset("genes", data=fixed(ids))
    g.create_dataset("gene_names", data=fixed(names))

# Two genome groups: not supported, must be a clear error.
with h5py.File(f"{out}/v2_two_genomes.h5", "w") as f:
    for genome in ("GRCh38", "mm10"):
        g = f.create_group(genome)
        g.create_dataset("barcodes", data=fixed(barcodes))
        g.create_dataset("data", data=data)
        g.create_dataset("indices", data=indices)
        g.create_dataset("indptr", data=indptr)
        g.create_dataset("shape", data=shape)
        g.create_dataset("genes", data=fixed(ids))
        g.create_dataset("gene_names", data=fixed(names))
