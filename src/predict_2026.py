"""Build a 2026 submission: 300 knockdowns x 400 cells x contexts A, B, C, raw counts.

Model ("transfer", tested in cross_line.py):
  1. For each target, average its share change (lfc, see cross_line.py) over every public
     cell line where it was switched off, on the genes that line measures. Genes no line
     measures get no change. Targets never switched off anywhere get the average change.
  2. Shrink the change by ALPHA (0.5 lowered error in 4 of 5 held-out lines in Exp 4).
     Optional denoising (K): first weight each line's change for gene g by z^2 / (z^2 + K),
     z = the change in counts / its noise (noise from src/noise.py), so changes that don't
     stand out from noise are pulled to zero.
     Optional per-knockdown weighting (KT): weight a line's whole change for a target by
     s / (s + KT), s = number of genes changed by more than 5 noise units (|z| > 5), not
     counting the switched-off gene itself. For the 2026 targets in K562 the median s is 1,
     so most of those profiles are close to pure noise.
  3. The switched-off gene itself always drops by the typical amount (not shrunk).
     Optional GAMMA: keep only GAMMA x the panel's average change (the part shared by all
     targets, mostly the source line's generic response) plus each target's own deviation
     from it. GAMMA = 1 changes nothing; 0 removes the shared part.
  4. Apply to each context's own control share and total counts.
  5. Draw 400 new cells per target: each cell's total is drawn from the context's control
     totals, then counts ~ gamma-Poisson around the predicted share, with per-gene
     overdispersion measured on the context's controls. No control cell is copied.

Arc's `vcc prep` loads the whole matrix (~33 GB for ~2e9 nonzeros), so this script
writes the same slim .h5ad in chunks and packages the .vcc (tar of meta.json +
pred.h5ad.zst) itself, matching vcc/prep.py `_write_vcc`.

Output: data/submissions/{name}.vcc
Usage: python src/predict_2026.py [name] [alpha] [gamma] [hc|none] [nb|-] [km]
       (default transfer_a05, ALPHA, 1, none, off, no map)
"""
import io
import json
import os
import sys
import tarfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zstandard as zstd

ROOT = Path(__file__).resolve().parent.parent
CTRL_DIR = Path(os.environ.get("VCC_CTRL_DIR", ROOT / "data/vcc/controls"))  # final round: set VCC_CTRL_DIR
SOURCES = ["h1", "k562", "rpe1", "hepg2", "jurkat"]
CONTEXTS = os.environ.get("VCC_CONTEXTS", "A,B,C").split(",")  # final round: e.g. VCC_CONTEXTS=D,E,F
CELLS = 400
ALPHA = 0.5
EPS = 1e-5
SEED = 0


def source_changes(genes, sources=SOURCES, k=None, kt=None):
    """Per target: summed lfc and number of lines measuring each gene, on the given genes.
    With k, each lfc is first weighted by z^2 / (z^2 + k); with kt, each target's whole
    lfc by s / (s + kt) (see module docstring)."""
    gidx = {g: i for i, g in enumerate(genes)}
    lfc_sum, lfc_n, ltr, own = {}, {}, {}, []
    line_means = []
    for name in sources:
        d = np.load(ROOT / "data/lines" / f"{name}.npz")
        keep = np.array([g in gidx for g in d["genes"]])
        cols = np.array([gidx[g] for g in d["genes"][keep]])
        counts, ctrl = d["counts"][:, keep], d["ctrl"][keep]
        share = counts / counts.sum(1, keepdims=True)
        ctrl_share = ctrl / ctrl.sum()
        lfc = np.log((share + EPS) / (ctrl_share + EPS))
        r = np.log(counts.sum(1) / ctrl.sum())
        local_all = {g: i for i, g in enumerate(d["genes"][keep])}
        if k is not None or kt is not None:
            nz = np.load(ROOT / "data/lines" / f"{name}_noise.npz")
            nidx = {g: i for i, g in enumerate(nz["genes"])}
            var = nz["var"][[nidx[g] for g in d["genes"][keep]]]
            diff = counts - ctrl[None] * np.exp(r)[:, None]
            z2 = diff ** 2 / (np.maximum(var, 1e-6)[None] * (1 / d["n"][:, None] + 1 / d["ctrl_n"]))
            if k is not None:
                lfc = lfc * z2 / (z2 + k)
            if kt is not None:
                z2 = z2.copy()
                for i, t in enumerate(d["targets"]):
                    if t in local_all:
                        z2[i, local_all[t]] = 0  # the switched-off gene itself doesn't count
                strength = (z2 > 25).sum(1)
                lfc = lfc * (strength / (strength + kt))[:, None]
        local = {g: i for i, g in enumerate(d["genes"][keep])}
        full_mean = np.zeros(len(genes))
        full_mean[cols] = lfc.mean(0)
        line_means.append((full_mean, cols, r.mean()))
        for i, t in enumerate(d["targets"]):
            if t not in lfc_sum:
                lfc_sum[t], lfc_n[t], ltr[t] = np.zeros(len(genes)), np.zeros(len(genes)), []
            lfc_sum[t][cols] += lfc[i]
            lfc_n[t][cols] += 1
            ltr[t].append(r[i])
            if t in local:
                own.append(lfc[i, local[t]])
    # average change over lines (each line weighted equally), for never-seen targets
    msum, mn = np.zeros(len(genes)), np.zeros(len(genes))
    for full_mean, cols, _ in line_means:
        msum[cols] += full_mean[cols]
        mn[cols] += 1
    mean_lfc = np.where(mn > 0, msum / np.maximum(mn, 1), 0)
    mean_ltr = np.mean([x for _, _, x in line_means])
    return lfc_sum, lfc_n, ltr, mean_lfc, mean_ltr, float(np.median(own))


def target_change(t, genes_idx, src, alpha=ALPHA):
    lfc_sum, lfc_n, ltr, mean_lfc, mean_ltr, own_drop = src
    if t in lfc_sum:
        lfc = np.where(lfc_n[t] > 0, lfc_sum[t] / np.maximum(lfc_n[t], 1), 0)
        r = np.mean(ltr[t])
    else:
        lfc, r = mean_lfc.copy(), mean_ltr
    lfc, r = alpha * lfc, alpha * r
    if t in genes_idx:
        lfc[genes_idx[t]] = own_drop
    return lfc, r


def recenter(changes, genes_idx, gamma=1.0):
    """Scale the panel's shared change by gamma: lfc_t -> lfc_t - (1 - gamma) * mean_t(lfc_t).
    The switched-off gene's own drop is left as it was."""
    if gamma == 1.0:
        return changes
    shared = np.mean([lfc for lfc, _ in changes.values()], axis=0)
    out = {}
    for t, (lfc, r) in changes.items():
        new = lfc - (1 - gamma) * shared
        if t in genes_idx:
            new[genes_idx[t]] = lfc[genes_idx[t]]
        out[t] = (new, r)
    return out


def line_common(genes, line, exclude=()):
    """A line's typical change: mean lfc (and log total ratio) over its knockdowns, skipping
    `exclude`, on the given genes (0 where the line doesn't measure a gene)."""
    d = np.load(ROOT / "data/lines" / f"{line}.npz")
    gidx = {g: i for i, g in enumerate(genes)}
    keep = np.array([g in gidx for g in d["genes"]])
    rows = ~np.isin(d["targets"], list(exclude))
    counts, ctrl = d["counts"][rows][:, keep], d["ctrl"][keep]
    share = counts / counts.sum(1, keepdims=True)
    lfc = np.log((share + EPS) / (ctrl / ctrl.sum() + EPS))
    out = np.zeros(len(genes))
    out[[gidx[g] for g in d["genes"][keep]]] = lfc.mean(0)
    return out, float(np.log(counts.sum(1) / ctrl.sum()).mean())


def add_common(changes, genes_idx, common, weight=1.0):
    """Add weight x a line's typical change to every target (own gene's drop kept) and use
    that line's typical total-count change."""
    lfc_c, r_c = common
    out = {}
    for t, (lfc, r) in changes.items():
        new = lfc + weight * lfc_c
        if t in genes_idx:
            new[genes_idx[t]] = lfc[genes_idx[t]]
        out[t] = (new, weight * r_c)
    return out


def k562_to_h1(genes, targets, lam=10.0, exclude=()):
    """Linear map (ridge) from a target's K562 change to its H1 change, trained on targets
    switched off in both lines (minus `exclude`). Returns {target: predicted H1-specific change}
    on `genes` (0 off the shared genes) for targets switched off in K562."""
    k, h = (np.load(ROOT / "data/lines" / f"{n}.npz") for n in ("k562", "h1"))
    common = np.intersect1d(np.intersect1d(k["genes"], h["genes"]), genes)

    def lfc(d):
        ix = {g: i for i, g in enumerate(d["genes"])}
        cols = np.array([ix[g] for g in common])
        c, ctrl = d["counts"][:, cols], d["ctrl"][cols]
        return np.log((c / c.sum(1, keepdims=True) + EPS) / (ctrl / ctrl.sum() + EPS))

    LK, LH = lfc(k), lfc(h)
    LK -= LK.mean(0)
    rk = {t: i for i, t in enumerate(k["targets"])}
    rh = {t: i for i, t in enumerate(h["targets"])}
    train = [t for t in h["targets"] if t in rk and t not in set(exclude)]
    X = LK[[rk[t] for t in train]]
    Y = LH[[rh[t] for t in train]]
    Y = Y - Y.mean(0)
    A = np.linalg.solve(X @ X.T + lam * np.eye(len(train)), Y)  # W = X^T A (dual form)
    gidx = {g: i for i, g in enumerate(genes)}
    cols = np.array([gidx[g] for g in common])
    out = {}
    for t in targets:
        if t in rk:
            v = np.zeros(len(genes))
            v[cols] = (LK[rk[t]] @ X.T) @ A
            out[t] = v
    print(f"k562_to_h1: trained on {len(train)} targets, {len(common):,} genes, {len(out)} predicted")
    return out


def add_mapped(changes, genes_idx, mapped):
    """Add the mapped change to each target's change (own gene's drop kept)."""
    out = {}
    for t, (lfc, r) in changes.items():
        new = lfc + mapped[t] if t in mapped else lfc
        if t in genes_idx:
            new[genes_idx[t]] = lfc[genes_idx[t]]
        out[t] = (new, r)
    return out


def add_neighbours(changes, genes, own_drop, near=1e3, far=5e3):
    """Switching a gene off also lowers genes that start right next to it on the DNA (H1: 0.28x
    within 1 kb, 0.40x within 5 kb; K562 shows the same, weaker). Set genes starting within
    `near` of the target's start to the typical own drop, and within `far` to half of it."""
    tss = pd.read_csv(ROOT / "data/annot/tss.csv", index_col=0).reindex(genes)
    chrom, pos = tss["chr"].to_numpy(), tss["tss"].to_numpy()
    out = {}
    for t, (lfc, r) in changes.items():
        new = lfc.copy()
        if t in tss.index and isinstance(tss.loc[t, "chr"], str):
            d = np.abs(pos - tss.loc[t, "tss"])
            same = (chrom == tss.loc[t, "chr"]) & (genes != t)
            new[same & (d < far)] = np.minimum(new[same & (d < far)], own_drop / 2)
            new[same & (d < near)] = np.minimum(new[same & (d < near)], own_drop)
        out[t] = (new, r)
    return out


def context_stats(X):
    """Control share, per-cell totals and per-gene overdispersion from control cells X."""
    X = sp.csr_matrix(X, dtype=np.float64)
    totals = np.asarray(X.sum(1)).ravel()
    share = np.asarray(X.sum(0)).ravel() / totals.sum()
    # per-gene overdispersion on depth-normalized counts: var = m + phi * m^2
    Y = sp.diags(np.median(totals) / totals) @ X
    m = np.asarray(Y.mean(0)).ravel()
    v = np.asarray(Y.multiply(Y).mean(0)).ravel() - m ** 2
    phi = np.clip((v - m) / np.maximum(m, 1e-12) ** 2, 0, 10)
    phi[m == 0] = 0
    return share, totals, phi


def sample_cells(share, totals, phi, ltr, rng, n=CELLS):
    lib = rng.choice(totals, n) * np.exp(ltr)
    mu = lib[:, None] * share[None, :]
    over = phi > 0
    g = np.ones_like(mu)
    shape = 1 / phi[over]
    g[:, over] = rng.gamma(shape, 1 / shape, size=(n, over.sum()))
    return sp.csr_matrix(rng.poisson(mu * g).astype(np.float32))


def write_h5ad(path, genes, targets, gen):
    """Slim h5ad like `vcc prep`: obs {target_gene, context}, var index genes, X float32 CSR."""
    n = len(CONTEXTS) * len(targets) * CELLS
    obs = pd.DataFrame({
        "target_gene": np.tile(np.repeat(targets, CELLS), len(CONTEXTS)),
        "context": np.repeat(CONTEXTS, len(targets) * CELLS),
    }, index=np.arange(n).astype(str))
    ad.AnnData(obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(path)
    nnz = 0
    with h5py.File(path, "a") as f:
        if "X" in f:
            del f["X"]
        first = True
        for block in gen:
            nnz += block.nnz
            if first:
                ad.io.write_elem(f, "X", block)
                X = ad.io.sparse_dataset(f["X"])
                first = False
            else:
                X.append(block)
    return n, nnz


def package(h5ad_path, vcc_path, n_obs, n_vars, nnz):
    zst = h5ad_path.with_suffix(".h5ad.zst")
    cctx = zstd.ZstdCompressor(level=3, threads=os.cpu_count() or 1)
    with open(h5ad_path, "rb") as src, open(zst, "wb") as dst:
        cctx.copy_stream(src, dst)
    meta = json.dumps({"cli_version": "0.2.2", "n_obs": n_obs, "n_vars": n_vars,
                       "nnz": nnz, "schema": 1}, sort_keys=True).encode()

    def norm(ti):
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        ti.mtime = 0
        return ti

    with tarfile.open(vcc_path, "w") as tar:
        ti = norm(tarfile.TarInfo("meta.json"))
        ti.size = len(meta)
        tar.addfile(ti, io.BytesIO(meta))
        tar.add(zst, arcname="pred.h5ad.zst", filter=norm)
    zst.unlink()


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "transfer_a05"
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else ALPHA
    gamma = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    hc = float(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != "none" else None
    nb = len(sys.argv) > 5 and sys.argv[5] == "nb"
    km = float(sys.argv[6]) if len(sys.argv) > 6 else None
    genes = pd.read_csv(CTRL_DIR / "gene_names.csv")["gene_name"].to_numpy(str)
    targets = pd.read_csv(CTRL_DIR / "pert_counts.csv")["target_gene"].to_numpy(str)
    gidx = {g: i for i, g in enumerate(genes)}
    src = source_changes(genes)
    seen = sum(t in src[0] for t in targets)
    print(f"targets switched off in a public line: {seen}/{len(targets)}; "
          f"typical own drop: {np.exp(src[5]):.2f}x")

    rng = np.random.default_rng(SEED)
    changes = recenter({t: target_change(t, gidx, src, alpha) for t in targets}, gidx, gamma)
    if hc is not None:  # H1's typical change (see MODEL.md)
        changes = add_common(changes, gidx, line_common(genes, "h1"), hc)
    if km is not None:  # K562 -> H1 linear map, ridge strength km (see k562_to_h1)
        changes = add_mapped(changes, gidx, k562_to_h1(genes, targets, km))
    if nb:  # neighbouring genes drop too (see add_neighbours)
        changes = add_neighbours(changes, genes, src[5])

    def blocks():
        for c in CONTEXTS:
            a = ad.read_h5ad(CTRL_DIR / f"context_{c}.h5ad")
            assert list(a.var_names) == list(genes)
            share, totals, phi = context_stats(a.X)
            for k, t in enumerate(targets):
                lfc, r = changes[t]
                s = np.clip((share + EPS) * np.exp(lfc) - EPS, 0, None)
                yield sample_cells(s / s.sum(), totals, phi, r, rng)
                if (k + 1) % 50 == 0:
                    print(f"  context {c}: {k + 1}/{len(targets)} targets", flush=True)

    out = ROOT / "data/submissions"
    out.mkdir(parents=True, exist_ok=True)
    h5 = out / f"{name}.h5ad"
    n_obs, nnz = write_h5ad(h5, genes, targets, blocks())
    print(f"wrote {h5} ({n_obs:,} cells, {nnz:,} nonzeros, {nnz / n_obs:,.0f} per cell)")
    package(h5, out / f"{name}.vcc", n_obs, len(genes), nnz)
    print(f"packaged {out / f'{name}.vcc'}")


if __name__ == "__main__":
    main()
