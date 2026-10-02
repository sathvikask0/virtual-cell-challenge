"""Which source transfers best to which cell line? Per-target correlation of each source's centered change (E_c)
with the real centered change in a held-out local-eval line (real.h5ad), over well-expressed genes.
Usage (from src/): python source_match.py LINE [LINE ...]   (needs data/atlas_shift_x from context_weights.py build)
"""
import sys
import anndata as ad
import numpy as np
import scipy.sparse as sp
import atlas_shift as A
import context_weights as CW

CW.use_x()
for line in sys.argv[1:]:
    a = ad.read_h5ad(A.ROOT / "data/local_eval" / line / "real.h5ad")
    lab = a.obs["target"].astype(str).to_numpy()
    genes = np.array(a.var_names, str)
    X = sp.csr_matrix(a.X)
    groups = sorted(set(lab))
    S = np.vstack([np.asarray(X[lab == g].sum(0)).ravel() for g in groups])
    c = S[groups.index("non-targeting")]; p0 = c / c.sum()
    targets = [g for g in groups if g != "non-targeting"]
    T = S[[groups.index(t) for t in targets]]
    p = (T + A.PRIOR * p0) / (T.sum(1, keepdims=True) + A.PRIOR)
    ec = np.log2((1e6 * p + 1) / (1e6 * p0 + 1))
    gi = {g: i for i, g in enumerate(genes)}
    for i, t in enumerate(targets):
        if t in gi: ec[i, gi[t]] = np.nan
    ec -= np.nanmean(ec, 0)
    expressed = 1e6 * p0 > 20
    srcs = dict(A.WEIGHTS); srcs["cd4"] = 1
    print(f"== {line}: {len(targets)} targets, {expressed.sum()} genes >20 CPM")
    for s in srcs:
        d = np.load(CW.C.PATH if s == "cd4" else A.OUT / f"{s}.npz", allow_pickle=True)
        sg = d["genes"].astype(str) if s == "cd4" else A.genes26()
        sgi = {g: i for i, g in enumerate(sg)}
        cols = np.array([sgi.get(g, -1) for g in genes])
        rows = {t: i for i, t in enumerate(d["targets"].astype(str))}
        r = []
        for i, t in enumerate(targets):
            if t not in rows: continue
            v = np.where(cols >= 0, d["ec"][rows[t]][np.maximum(cols, 0)], np.nan)
            m = expressed & np.isfinite(v) & np.isfinite(ec[i])
            if m.sum() > 100: r.append(np.corrcoef(v[m], ec[i, m])[0, 1])
        print(f"  {s:8s} n={len(r):3d}  median r={np.median(r):.3f}  mean r={np.mean(r):.3f}")
