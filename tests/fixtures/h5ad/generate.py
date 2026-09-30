#!/usr/bin/env python3
"""Regenerates the h5ad fixtures used by tests/h5ad.rs.

Every fixture encodes the same 2 cells x 3 genes matrix unless noted:

    cells      ENSG1  ENSG2  ENSG3
    AAAC-1       5      0      1
    AAAG-1       0      3      0

Requires: anndata >= 0.10, h5py, numpy, scipy.

    python3 tests/fixtures/h5ad/generate.py tests/fixtures/h5ad
"""
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

out = sys.argv[1] if len(sys.argv) > 1 else "."
X = np.array([[5, 0, 1], [0, 3, 0]], dtype=np.float32)
obs = pd.DataFrame(index=["AAAC-1", "AAAG-1"])
var = pd.DataFrame({"gene_symbols": ["A", "B", "C"]}, index=["ENSG1", "ENSG2", "ENSG3"])


def write(name, adata):
    adata.write_h5ad(f"{out}/{name}.h5ad")


# Modern anndata layouts: index and string columns are nullable-string-array groups.
write("csr", ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var))
write("csc", ad.AnnData(X=sp.csc_matrix(X), obs=obs, var=var))
write("dense", ad.AnnData(X=X, obs=obs, var=var))

# Symbols only as a categorical column (CELLxGENE-style feature_name); int64 data.
var_cat = pd.DataFrame({"feature_name": pd.Categorical(["A", "B", "C"])}, index=var.index)
write("categorical_int64", ad.AnnData(X=sp.csr_matrix(X.astype(np.int64)), obs=obs, var=var_cat))

# Custom index column names declared through the `_index` attribute.
a = ad.AnnData(X=sp.csr_matrix(X), obs=obs.copy(), var=var.copy())
a.obs.index.name = "barcode"
a.var.index.name = "gene_ids"
write("named_index", a)

# Normalized X with raw counts kept under raw/X.
a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
a.raw = a.copy()
a.X = sp.csr_matrix(np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4))
write("raw_layer", a)

# CSR with an explicit zero, unsorted indices and a duplicate coordinate:
#   row 0: (gene2, 1) (gene0, 5) (gene2, 0) (gene0, 1)  -> gene0 = 6, gene2 = 1
#   row 1: (gene1, 3)
write("messy", ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var))
with h5py.File(f"{out}/messy.h5ad", "r+") as f:
    for k in ("data", "indices", "indptr"):
        del f[f"X/{k}"]
    f["X"].create_dataset("data", data=np.array([1, 5, 0, 1, 3], dtype=np.float32))
    f["X"].create_dataset("indices", data=np.array([2, 0, 2, 0, 1], dtype=np.int32))
    f["X"].create_dataset("indptr", data=np.array([0, 4, 5], dtype=np.int32))

# NaN in X/data.
write("nonfinite", ad.AnnData(
    X=sp.csr_matrix(np.array([[np.nan, 1, 0], [0, 2, 0]], dtype=np.float32)), obs=obs, var=var))

# First stored index points past the last gene.
write("out_of_range", ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var))
with h5py.File(f"{out}/out_of_range.h5ad", "r+") as f:
    idx = f["X/indices"][...]
    idx[0] = 7
    del f["X/indices"]
    f["X"].create_dataset("indices", data=idx)

# indptr one entry short: structural corruption.
write("short_indptr", ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var))
with h5py.File(f"{out}/short_indptr.h5ad", "r+") as f:
    ip = f["X/indptr"][...]
    del f["X/indptr"]
    f["X"].create_dataset("indptr", data=ip[:-1])


def write_legacy(path, fixed_width):
    """anndata <= 0.10 layout: obs/var groups whose index is a plain string dataset."""
    with h5py.File(path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array([2, 3], dtype=np.int64)
        m = sp.csr_matrix(X)
        g.create_dataset("data", data=m.data.astype(np.float32))
        g.create_dataset("indices", data=m.indices.astype(np.int32))
        g.create_dataset("indptr", data=m.indptr.astype(np.int32))
        for grp, names in (("obs", ["AAAC-1", "AAAG-1"]), ("var", ["ENSG1", "ENSG2", "ENSG3"])):
            gg = f.create_group(grp)
            gg.attrs["_index"] = "_index"
            gg.attrs["encoding-type"] = "dataframe"
            gg.attrs["encoding-version"] = "0.1.0"
            gg.attrs["column-order"] = np.array([], dtype="S1")
            dtype = "S12" if fixed_width else h5py.string_dtype()
            gg.create_dataset("_index", data=np.array(names, dtype=dtype))
        f["var"].create_dataset("gene_symbols", data=np.array(["A", "B", "C"], dtype="S4"))


write_legacy(f"{out}/legacy_vlen.h5ad", fixed_width=False)
write_legacy(f"{out}/legacy_fixed.h5ad", fixed_width=True)

# 10x Feature Barcoding as written by scanpy.read_10x_h5: var/feature_types is
# categorical, one antibody feature with counts far above the gene scale.
var_fb = pd.DataFrame(
    {
        "gene_symbols": ["A", "B", "CD3_TotalSeqB"],
        "feature_types": pd.Categorical(["Gene Expression", "Gene Expression", "Antibody Capture"]),
    },
    index=["ENSG1", "ENSG2", "CD3_TotalSeqB"],
)
X_fb = np.array([[5, 0, 9000], [0, 3, 8000]], dtype=np.float32)
write("feature_types", ad.AnnData(X=sp.csr_matrix(X_fb), obs=obs, var=var_fb))
