"""K562 genome-wide CRISPRi from the per-cell file, into the atlas_shift cache (data/atlas_shift/k562_cpm.npz).

Why: our `k562` source uses the averaged (bulk) file, so the per-cell change E_c is approximated from pooled
shares. Here both are computed properly, as in the public #82 prepare_k562 (kaipengm2, MIT):
  - per target: summed counts, and the mean of per-cell normalized counts (CPM)
  - controls per sequencing batch (gem_group); each target is compared with controls from the batches its
    cells came from (weighted by cells for CPM, by library size for pooled shares)
  - a 1e5-count prior toward the matched control, as in atlas_shift
Then E_c = log2((mean CPM_t + 1) / (matched control CPM + 1)), E_b = log1p(5e4 p_t) - log1p(5e4 p_0),
centered over all K562 targets (own gene excluded), and stored for the 2026 + H1 targets on the 2026 gene list.

Input: data/replogle/K562_gwps_raw_singlecell_01.h5ad (66 GB, figshare 35775507).
Usage: python src/k562_cells.py
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_shift as A

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data/replogle/K562_gwps_raw_singlecell_01.h5ad"
CTRL = "non-targeting"
CHUNK = 4096


def main():
    a = ad.read_h5ad(SRC, backed="r")
    labels = a.obs["gene"].astype(str).to_numpy()
    bcode, batches = pd.factorize(a.obs["gem_group"].astype(str))
    genes0 = a.var["gene_name"].astype(str).to_numpy()
    genes = np.array(list(dict.fromkeys(genes0)), dtype=str)
    gl = {g: i for i, g in enumerate(genes)}
    proj = sp.csr_matrix((np.ones(len(genes0)), (np.arange(len(genes0)), [gl[g] for g in genes0])),
                         shape=(len(genes0), len(genes)))
    targets = np.array(sorted(set(labels) - {CTRL}))
    tl = {t: i for i, t in enumerate(targets)}
    T, B = len(targets), len(batches)
    group = np.array([tl.get(t, -1) for t in labels])
    ctrl = labels == CTRL
    group[ctrl] = T + bcode[ctrl]
    rows = np.flatnonzero(group >= 0)
    counts = np.zeros((T + B, len(genes)))
    cpm = np.zeros_like(counts)
    n = np.bincount(group[rows], minlength=T + B)
    tb_n = np.zeros((T, B))
    tb_lib = np.zeros((T, B))
    for k, lo in enumerate(range(0, len(rows), CHUNK)):
        sel = rows[lo:lo + CHUNK]
        x = sp.csr_matrix(a.X[sel]).astype(np.float64) @ proj
        depth = np.asarray(x.sum(1)).ravel()
        gi = group[sel]
        asg = sp.csr_matrix((np.ones(len(sel)), (gi, np.arange(len(sel)))), shape=(T + B, len(sel)))
        counts += (asg @ x).toarray()
        cpm += (asg @ sp.diags(1e6 / depth) @ x).toarray()
        p = gi < T
        np.add.at(tb_n, (gi[p], bcode[sel][p]), 1)
        np.add.at(tb_lib, (gi[p], bcode[sel][p]), depth[p])
        if k % 100 == 0:
            print(f"  {lo + len(sel):,}/{len(rows):,} cells", flush=True)
    a.file.close()

    cprob = counts[T:] / counts[T:].sum(1, keepdims=True)          # per-batch control shares
    ccpm = cpm[T:] / n[T:, None]                                   # per-batch control mean CPM
    nt = n[:T]
    w_cells = tb_n / np.maximum(nt[:, None], 1)
    w_lib = tb_lib / np.maximum(tb_lib.sum(1, keepdims=True), 1e-9)
    p0 = w_lib @ cprob                                             # matched control shares, per target
    c0 = w_cells @ ccpm                                            # matched control mean CPM, per target
    sums = counts[:T]
    p = (sums + A.PRIOR * p0) / (sums.sum(1, keepdims=True) + A.PRIOR)
    frac = sums.sum(1) / (sums.sum(1) + A.PRIOR)
    mcpm = frac[:, None] * cpm[:T] / np.maximum(nt[:, None], 1) + (1 - frac[:, None]) * c0
    ec = np.log2((mcpm + 1) / (c0 + 1))
    eb = np.log1p(5e4 * p) - np.log1p(5e4 * p0)
    ok = nt >= 20
    own = np.zeros_like(ec, bool)
    for i, t in enumerate(targets):
        if t in gl:
            own[i, gl[t]] = True
    m = ok[:, None] & ~own
    ec -= np.where(m, ec, 0).sum(0) / m.sum(0)
    eb -= np.where(m, eb, 0).sum(0) / m.sum(0)

    g26 = A.genes26()
    gidx = {g: i for i, g in enumerate(g26)}
    cols = np.array([gidx.get(g, -1) for g in genes])
    want = A.wanted_targets()
    keep = np.array([ok[i] and t in want for i, t in enumerate(targets)])
    EC = np.full((keep.sum(), len(g26)), np.nan, np.float32)
    EB = EC.copy()
    EC[:, cols[cols >= 0]] = ec[keep][:, cols >= 0]
    EB[:, cols[cols >= 0]] = eb[keep][:, cols >= 0]
    np.savez(A.OUT / "k562_cpm.npz", targets=targets[keep], ec=EC, eb=EB)
    print(f"k562_cpm: {T:,} targets, {keep.sum()} cached (2026 + H1), {len(genes):,} genes, "
          f"{B} batches, median {np.median(nt[ok]):.0f} cells per target", flush=True)


if __name__ == "__main__":
    main()
