"""Per-gene noise of each line: variance of raw counts per control cell.

Used to tell a real knockdown effect from noise: the standard error of a target's average
count for gene g is about sqrt(var[g] / n_cells).

  per-cell files (h1, hepg2, jurkat): variance across control cells
  Replogle averaged files (k562, rpe1): spread of the ~100-600 control-guide averages,
      var = sum_j n_j (mean_j - overall)^2 / (J - 1), since each guide average has variance var / n_j.
      This also picks up guide-to-guide and batch noise, which is realistic.

Output: data/lines/{line}_noise.npz with genes, var (per control cell).
Usage: python src/noise.py [line ...]
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
CTRL = "non-targeting"
CHUNK = 20_000
PER_CELL = {
    "h1": (ROOT / "data/vcc/train/adata_Training.h5ad", "target_gene"),
    "hepg2": (ROOT / "data/nadig/hepg2_raw_singlecell.h5ad", "gene"),
    "jurkat": (ROOT / "data/nadig/jurkat_raw_singlecell.h5ad", "gene"),
}
BULK = {
    "k562": ROOT / "data/replogle/K562_gwps_raw_bulk_01.h5ad",
    "rpe1": ROOT / "data/replogle/rpe1_raw_bulk_01.h5ad",
}


def per_cell(path, col):
    a = ad.read_h5ad(path, backed="r")
    is_ctrl = (a.obs[col].astype(str) == CTRL).to_numpy()
    s = ss = 0
    n = 0
    for start in range(0, a.n_obs, CHUNK):
        keep = is_ctrl[start:start + CHUNK]
        if not keep.any():
            continue
        X = sp.csr_matrix(a.X[start:start + CHUNK], dtype=np.float64)[keep]
        s = s + np.asarray(X.sum(0)).ravel()
        ss = ss + np.asarray(X.multiply(X).sum(0)).ravel()
        n += keep.sum()
    mean = s / n
    genes = np.array(a.var["gene_name"] if "gene_name" in a.var else a.var_names, dtype=str)
    return genes, ss / n - mean ** 2


def bulk(path):
    a = ad.read_h5ad(path)
    c = np.array(["non-targeting" in i for i in a.obs_names])
    w = a.obs["num_cells_filtered"].fillna(a.obs["num_cells_unfiltered"]).to_numpy(float)[c]
    X = np.asarray(a.X, dtype=float)[c]
    mean = (X * w[:, None]).sum(0) / w.sum()
    var = (w[:, None] * (X - mean) ** 2).sum(0) / (len(w) - 1)
    return np.array(a.var["gene_name"], dtype=str), var


def main():
    for line in sys.argv[1:] or [*PER_CELL, *BULK]:
        genes, var = per_cell(*PER_CELL[line]) if line in PER_CELL else bulk(BULK[line])
        _, first = np.unique(genes, return_index=True)
        first = np.sort(first)
        np.savez_compressed(ROOT / "data/lines" / f"{line}_noise.npz", genes=genes[first],
                            var=var[first])
        print(f"{line:7} median per-cell variance {np.median(var):.3f}")


if __name__ == "__main__":
    main()
