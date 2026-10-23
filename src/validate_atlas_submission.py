"""Bounded-memory validation of a full generated h5ad and its custom .vcc package."""
import hashlib
import json
import sys
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zstandard as zstd
from anndata.io import read_elem

from atlas_shift import ROOT, ctrl_dir, genes26


def digest(stream):
    sha = hashlib.sha256()
    for block in iter(lambda: stream.read(8 << 20), b""):
        sha.update(block)
    return sha.hexdigest()


def main(name):
    out = ROOT / "data/submissions"
    h5, archive = out / f"{name}.h5ad", out / f"{name}.vcc"
    import predict_2026 as P
    expected_targets = set(pd.read_csv(ctrl_dir() / "pert_counts.csv").target_gene)
    genes = genes26()
    with h5py.File(h5) as f:
        obs, var = read_elem(f["obs"]), read_elem(f["var"])
        if list(var.index) != list(genes) or not obs.index.is_unique:
            raise ValueError("Gene axis or cell IDs invalid")
        if set(obs.context.astype(str)) != set(P.CONTEXTS):
            raise ValueError("Context panel differs")
        for context, group in obs.groupby("context", observed=True):
            counts = group.target_gene.astype(str).value_counts()
            if set(counts.index) != expected_targets or not np.all(counts.to_numpy() == 400):
                raise ValueError(f"Incomplete target/cell panel: {context}")
        x = f["X"]
        if x.attrs["encoding-type"] != "csr_matrix" or tuple(x.attrs["shape"]) != (360000, len(genes)):
            raise ValueError("Sparse matrix schema differs")
        ptr = x["indptr"][:]
        nnz = len(x["data"])
        if len(x["indices"]) != nnz or len(ptr) != 360001 or ptr[0] != 0 or ptr[-1] != nnz or np.any(np.diff(ptr) < 0):
            raise ValueError("CSR pointers invalid")
        for start in range(0, nnz, 2_000_000):
            values = x["data"][start:start+2_000_000]
            indices = x["indices"][start:start+2_000_000]
            if not np.isfinite(values).all() or np.any(values < 0) or np.any(values != np.floor(values)):
                raise ValueError("Counts must be finite nonnegative integers")
            if np.any(indices < 0) or np.any(indices >= len(genes)):
                raise ValueError("CSR gene index out of range")
    print(f"Validated full panel, CSR structure and all {nnz:,} stored counts", flush=True)
    with h5.open("rb") as stream:
        h5_sha = digest(stream)
    with tarfile.open(archive) as tar:
        if sorted(tar.getnames()) != ["meta.json", "pred.h5ad.zst"]:
            raise ValueError("Unexpected archive members")
        meta = json.load(tar.extractfile("meta.json"))
        if meta.get("schema") != 1 or (meta.get("n_obs"),meta.get("n_vars"),meta.get("nnz")) != (360000,len(genes),nnz):
            raise ValueError("Archive metadata differs from h5ad")
        with zstd.ZstdDecompressor().stream_reader(tar.extractfile("pred.h5ad.zst")) as stream:
            packed_sha = digest(stream)
        if packed_sha != h5_sha:
            raise ValueError("Packed h5ad checksum differs")
    with archive.open("rb") as stream:
        archive_sha = digest(stream)
    report = dict(name=name, n_obs=360000, n_vars=len(genes), nnz=nnz,
                  targets_per_context=300, cells_per_target=400, contexts=list(P.CONTEXTS),
                  h5ad_sha256=h5_sha, vcc_sha256=archive_sha, checks="full counts, axes, CSR, panel, archive metadata and decompressed payload hash",
                  limitation="Local format/integrity validation; no leaderboard upload or score")
    (out / f"validation_{name}.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
