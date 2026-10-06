"""Put every cell line into one format: average raw counts per cell, per switched-off gene.

Lines with answers (used to practice predicting an unseen cell line):
  h1      VCC 2025 H1 stem cells       data/pseudobulk/{train,val,test}.npz  (300 targets)
  k562    Replogle 2022 genome-wide    data/replogle/K562_gwps_raw_bulk_01.h5ad  (~9.9k targets)
  rpe1    Replogle 2022 essential      data/replogle/rpe1_raw_bulk_01.h5ad  (~2.4k targets)
  hepg2   Nadig 2025 essential         data/nadig/hepg2_raw_singlecell.h5ad
  jurkat  Nadig 2025 essential         data/nadig/jurkat_raw_singlecell.h5ad
Lines without answers (2026 challenge, controls only):
  ctx_A, ctx_B, ctx_C                  data/vcc/controls/context_{A,B,C}.h5ad

Replogle bulk files only keep ~8k well-measured genes, so totals there cover those genes only.
Compare lines on their shared genes (see cross_line.py).

Output: data/lines/{line}.npz with
  genes   (G,)    gene symbols
  targets (T,)    switched-off gene symbols
  counts  (T, G)  average raw counts per cell
  n       (T,)    cells per target
  ctrl    (G,)    average raw counts per control cell
  ctrl_n  scalar  control cells

Usage: python src/lines.py [line ...]
"""
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
CTRL = "non-targeting"
CHUNK = 20_000


def dedup(genes, *arrays):
    """Keep the first column of each gene symbol."""
    _, first = np.unique(genes, return_index=True)
    first = np.sort(first)
    return (genes[first],) + tuple(a[..., first] for a in arrays)


def weighted_by_target(targets, X, w):
    """Average rows of X per target, weighted by cell counts w."""
    df = pd.DataFrame({"t": targets, "w": w})
    out_t, out_x, out_n = [], [], []
    for t, idx in df.groupby("t").indices.items():
        ww = w[idx]
        out_t.append(t)
        out_x.append((X[idx] * ww[:, None]).sum(0) / ww.sum())
        out_n.append(ww.sum())
    return np.array(out_t), np.array(out_x), np.array(out_n)


def from_h1():
    parts = [np.load(ROOT / "data/pseudobulk" / f"{s}.npz") for s in ["train", "val", "test"]]
    genes = parts[0]["genes"]
    t = np.concatenate([p["targets"] for p in parts])
    s = np.concatenate([p["sum_counts"] for p in parts])
    n = np.concatenate([p["n"] for p in parts]).astype(float)
    c = np.where(t == CTRL)[0][0]  # same control cells in every split file
    k = t != CTRL
    return genes, t[k], s[k] / n[k, None], n[k], s[c] / n[c], n[c]


def from_replogle_bulk(path):
    a = ad.read_h5ad(path)
    t = np.array([i.split("_")[1] for i in a.obs_names])
    # some rows have no filtered cell count; fall back to the unfiltered one
    w = a.obs["num_cells_filtered"].fillna(a.obs["num_cells_unfiltered"]).to_numpy(float)
    X = np.asarray(a.X, dtype=float)  # already average raw counts per cell
    targets, counts, n = weighted_by_target(t, X, w)
    k = targets != CTRL
    c = np.where(~k)[0][0]
    return (np.array(a.var["gene_name"], dtype=str), targets[k], counts[k], n[k],
            counts[c], n[c])


def from_singlecell(path, col):
    a = ad.read_h5ad(path, backed="r")
    cat = a.obs[col].astype(str).astype("category")
    codes, labels = cat.cat.codes.to_numpy(), np.array(cat.cat.categories, dtype=str)
    sums = np.zeros((len(labels), a.n_vars))
    for start in range(0, a.n_obs, CHUNK):
        stop = min(start + CHUNK, a.n_obs)
        X = sp.csr_matrix(a.X[start:stop], dtype=np.float64)
        G = sp.csr_matrix((np.ones(stop - start), (codes[start:stop], np.arange(stop - start))),
                          shape=(len(labels), stop - start))
        sums += np.asarray((G @ X).todense())
        print(f"  {path.name}: {stop:,}/{a.n_obs:,} cells", flush=True)
    n = np.bincount(codes, minlength=len(labels)).astype(float)
    means = sums / n[:, None]
    var = a.var
    genes = np.array(var["gene_name"] if "gene_name" in var else a.var_names, dtype=str)
    k = labels != CTRL
    c = np.where(~k)[0][0]
    return genes, labels[k], means[k], n[k], means[c], n[c]


LINES = {
    "h1": from_h1,
    "k562": lambda: from_replogle_bulk(ROOT / "data/replogle/K562_gwps_raw_bulk_01.h5ad"),
    "rpe1": lambda: from_replogle_bulk(ROOT / "data/replogle/rpe1_raw_bulk_01.h5ad"),
    "hepg2": lambda: from_singlecell(ROOT / "data/nadig/hepg2_raw_singlecell.h5ad", "gene"),
    "jurkat": lambda: from_singlecell(ROOT / "data/nadig/jurkat_raw_singlecell.h5ad", "gene"),
    **{f"ctx_{c}": (lambda c=c: from_singlecell(
        ROOT / f"data/vcc/controls/context_{c}.h5ad", "target_gene")) for c in "ABC"},
    **{f"ctx_{c}": (lambda c=c: from_singlecell(  # final round: VCC_CTRL_DIR + VCC_CONTEXTS
        Path(os.environ["VCC_CTRL_DIR"]) / f"context_{c}.h5ad", "target_gene"))
       for c in os.environ.get("VCC_CONTEXTS", "").split(",") if c and c not in "ABC"},
}


def main():
    out_dir = ROOT / "data/lines"
    out_dir.mkdir(parents=True, exist_ok=True)
    for line in sys.argv[1:] or LINES:
        genes, targets, counts, n, ctrl, ctrl_n = LINES[line]()
        genes, counts, ctrl = dedup(genes, counts, ctrl)
        np.savez_compressed(out_dir / f"{line}.npz", genes=genes, targets=targets,
                            counts=counts, n=n, ctrl=ctrl, ctrl_n=ctrl_n)
        print(f"{line:7} {len(targets):6,} targets  {len(genes):6,} genes  "
              f"control total {ctrl.sum():8,.0f} counts/cell  ({ctrl_n:,.0f} control cells)")


if __name__ == "__main__":
    main()
