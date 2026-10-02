"""KOLF2.1J iPSC CRISPRi atlas (Nature Biotech 2026, figshare 27261219, CC BY) as an atlas_shift source.

Input: data/kolf/KOLF_Strong_Perturbations.h5ad (47 GB; 1,656 targets with a significant response, 232k cells,
3 batches; raw counts in layers/counts, stored gene-major (CSC)).
Same statistics as src/k562_cells.py: per target, summed counts and mean per-cell CPM; controls (NTC) per batch,
matched to each target's batches; 1e5-count prior; centered over all KOLF targets with the own gene excluded.
Caveat: only strong knockdowns are in this file, so the centering average is over strong responses.

Writes kolf.npz into data/atlas_shift/ and data/atlas_shift_x/ (2026 + H1 + Jurkat/HepG2 local targets).
Usage: python src/kolf_cells.py
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

import atlas_shift as A
import context_weights as CW

SRC = A.ROOT / "data/kolf/KOLF_Strong_Perturbations.h5ad"
CTRL = "NTC"
GCHUNK = 2000


def main():
    f = h5py.File(SRC, "r")
    obs = read_elem(f["obs"])
    var_names = read_elem(f["var"]).index.astype(str).to_numpy()
    labels = obs["gene_target"].astype(str).to_numpy()
    bcode, batches = pd.factorize(obs["batch"].astype(str))
    targets = np.array(sorted(set(labels) - {CTRL}))
    tl = {t: i for i, t in enumerate(targets)}
    T, B = len(targets), len(batches)
    group = np.array([tl.get(t, -1) for t in labels])
    ctrl = labels == CTRL
    group[ctrl] = T + bcode[ctrl]
    n_cells = len(labels)
    asg = sp.csr_matrix((np.ones(n_cells), (group, np.arange(n_cells))), shape=(T + B, n_cells))

    c = f["layers/counts"]
    indptr = c["indptr"][:]
    G = len(indptr) - 1
    depth = np.zeros(n_cells)
    for lo in range(0, G, GCHUNK):  # pass 1: per-cell depth
        hi = min(lo + GCHUNK, G)
        a, b = indptr[lo], indptr[hi]
        depth += np.bincount(c["indices"][a:b], weights=c["data"][a:b], minlength=n_cells)
    inv = 1e6 / np.maximum(depth, 1)
    counts = np.zeros((T + B, G))
    cpm = np.zeros((T + B, G))
    for lo in range(0, G, GCHUNK):  # pass 2: group sums of counts and of per-cell CPM
        hi = min(lo + GCHUNK, G)
        a, b = indptr[lo], indptr[hi]
        X = sp.csc_matrix((c["data"][a:b].astype(np.float64), c["indices"][a:b], indptr[lo:hi + 1] - a),
                          shape=(n_cells, hi - lo))
        counts[:, lo:hi] = (asg @ X).toarray()
        cpm[:, lo:hi] = (asg @ sp.diags(inv) @ X).toarray()
        print(f"  genes {hi:,}/{G:,}", flush=True)
    f.close()

    # duplicate symbols: sum counts / cpm
    genes, inv_g = np.unique(var_names, return_inverse=True)
    proj = sp.csr_matrix((np.ones(G), (np.arange(G), inv_g)), shape=(G, len(genes)))
    counts, cpm = counts @ proj, cpm @ proj

    n = np.bincount(group, minlength=T + B)
    tb_n = np.zeros((T, B))
    tb_lib = np.zeros((T, B))
    p = group < T
    np.add.at(tb_n, (group[p], bcode[p]), 1)
    np.add.at(tb_lib, (group[p], bcode[p]), depth[p])
    cprob = counts[T:] / counts[T:].sum(1, keepdims=True)
    ccpm = cpm[T:] / n[T:, None]
    nt = n[:T]
    p0 = (tb_lib / np.maximum(tb_lib.sum(1, keepdims=True), 1e-9)) @ cprob
    c0 = (tb_n / np.maximum(nt[:, None], 1)) @ ccpm
    sums = counts[:T]
    pt = (sums + A.PRIOR * p0) / (sums.sum(1, keepdims=True) + A.PRIOR)
    frac = sums.sum(1) / (sums.sum(1) + A.PRIOR)
    mcpm = frac[:, None] * cpm[:T] / np.maximum(nt[:, None], 1) + (1 - frac[:, None]) * c0
    ec = np.log2((mcpm + 1) / (c0 + 1))
    eb = np.log1p(5e4 * pt) - np.log1p(5e4 * p0)
    gl = {g: i for i, g in enumerate(genes)}
    ok = nt >= 20
    m = ok[:, None] & np.ones_like(ec, bool)
    for i, t in enumerate(targets):
        if t in gl:
            m[i, gl[t]] = False
    ec -= np.where(m, ec, 0).sum(0) / m.sum(0)
    eb -= np.where(m, eb, 0).sum(0) / m.sum(0)

    g26 = A.genes26()
    gidx = {g: i for i, g in enumerate(g26)}
    cols = np.array([gidx.get(g, -1) for g in genes])
    want = A.wanted_targets() | CW.local_targets()
    keep = np.array([ok[i] and t in want for i, t in enumerate(targets)])
    EC = np.full((keep.sum(), len(g26)), np.nan, np.float32)
    EB = EC.copy()
    EC[:, cols[cols >= 0]] = ec[keep][:, cols >= 0]
    EB[:, cols[cols >= 0]] = eb[keep][:, cols >= 0]
    for out in (A.OUT, CW.XOUT):
        np.savez(out / "kolf.npz", targets=targets[keep], ec=EC, eb=EB)
    print(f"kolf: {T:,} targets, {keep.sum()} cached, {len(genes):,} genes, {B} batches, "
          f"median {np.median(nt[ok]):.0f} cells per target", flush=True)


if __name__ == "__main__":
    main()
