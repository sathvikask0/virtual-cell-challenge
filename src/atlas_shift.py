"""The public #82 recipe ("AtlasShift", kaipengm2, MIT), on our sources, testable on H1.

What differs from src/predict_2026.py:
  - two targets per knockdown instead of one:
      per-cell-normalized mean: m * 2^(AC * clip(E_c, 3)),  E_c = log2((1e6 p_t + 1) / (1e6 p_0 + 1))
      pooled profile:           log1p(5e4 q_t) = log1p(5e4 q) + AB * clip(E_b, 3),
                                E_b = log1p(5e4 p_t) - log1p(5e4 p_0)
    (the per-cell mean drives the significance tests, the pooled profile drives mse and pds,
    so they can be scaled separately: AC = 0.6, AB = 0.3 in #82)
  - each source's change is centered (its average over all its knockdowns removed, own gene excluded),
    then sources are averaged with weights (K562 2, H1 2, HCT116 1, HEK293T 1)
  - p_t gets a count prior: (sum_t + PRIOR * p_0) / (total_t + PRIOR)
  - neighbor prior: genes starting within 500 bp of the target keep 15%, ramping back to 100% at 5 kb
  - cells: 400 templates, each the average of POOL real control cells, reweighted so the group's
    per-cell mean and pooled profile hit both targets exactly, then rounded keeping each cell's total
    (third_party/atlasshift_model.dual_moment_counts). No random sampling.

Here K562 and X-Atlas only have per-target mean counts, so E_c uses pooled shares as a stand-in for
mean CPM (#82 uses per-cell CPM statistics).

Commands:
  python src/atlas_shift.py build                 cache centered source changes (data/atlas_shift/)
  python src/atlas_shift.py local LINE NAME [ac=0.6] [ab=0.3] [pool=4]   predict + score (as local_eval)
  python src/atlas_shift.py submission NAME [ac=0.6] [ab=0.3] [pool=4]   build data/submissions/NAME.vcc
"""
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from third_party.atlasshift_model import apply_promoter_prior, desired_mean, dual_moment_counts

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(os.environ.get("ATLAS_CACHE", ROOT / "data/atlas_shift"))  # final round: separate cache dir
WEIGHTS = {"k562": 2.0, "h1": 2.0, "hct116": 1.0, "hek293t": 1.0}
if os.environ.get("ATLAS_K562") == "cpm":  # K562 from the per-cell file (src/k562_cells.py)
    WEIGHTS = {("k562_cpm" if k == "k562" else k): w for k, w in WEIGHTS.items()}
PRIOR = 1e5
CTX_PROFILES = {}  # context -> profiles function, overrides `profiles` for that 2026 context
EXPR_W = None  # (dest cpm over genes_out, gamma): per-gene source weight min(1, 4*(src+1)/(dest+1))**gamma
CHUNK = 1000


def ctrl_dir():
    return Path(os.environ.get("VCC_CTRL_DIR", ROOT / "data/vcc/controls"))


def genes26():
    return pd.read_csv(ctrl_dir() / "gene_names.csv")["gene_name"].to_numpy(str)


def wanted_targets():
    t26 = set(pd.read_csv(ctrl_dir() / "pert_counts.csv")["target_gene"])
    return t26 | set(np.load(ROOT / "data/lines/h1.npz")["targets"])


def build():
    """Centered changes in both spaces for the 2026 + H1 targets, on the 2026 gene list."""
    OUT.mkdir(exist_ok=True)
    genes = genes26()
    gidx = {g: i for i, g in enumerate(genes)}
    want = wanted_targets()
    for line in WEIGHTS:
        d = np.load(ROOT / "data/lines" / f"{line}.npz")
        cols = np.array([gidx.get(g, -1) for g in d["genes"]])
        ok = cols >= 0
        ctrl = d["ctrl"][ok]
        p0 = ctrl / ctrl.sum()
        own = {t: np.where(d["genes"][ok] == t)[0] for t in d["targets"]}

        def effects(rows):
            sums = d["counts"][rows][:, ok] * d["n"][rows, None]
            p = (sums + PRIOR * p0) / (sums.sum(1, keepdims=True) + PRIOR)
            ec = np.log2((1e6 * p + 1) / (1e6 * p0 + 1))
            eb = np.log1p(5e4 * p) - np.log1p(5e4 * p0)
            return ec, eb

        T = len(d["targets"])
        sc, sb, cnt = np.zeros(ok.sum()), np.zeros(ok.sum()), np.zeros(ok.sum())
        for a in range(0, T, CHUNK):
            rows = np.arange(a, min(a + CHUNK, T))
            ec, eb = effects(rows)
            m = np.ones_like(ec, bool)
            for k, r in enumerate(rows):
                m[k, own[d["targets"][r]]] = False
            sc += np.where(m, ec, 0).sum(0)
            sb += np.where(m, eb, 0).sum(0)
            cnt += m.sum(0)
        rows = np.array([i for i, t in enumerate(d["targets"]) if t in want])
        ec, eb = effects(rows)
        ec -= sc / cnt
        eb -= sb / cnt
        EC = np.full((len(rows), len(genes)), np.nan, np.float32)
        EB = EC.copy()
        EC[:, cols[ok]] = ec
        EB[:, cols[ok]] = eb
        np.savez(OUT / f"{line}.npz", targets=d["targets"][rows], ec=EC, eb=EB)
        print(f"{line}: {len(rows)} targets cached", flush=True)


def fused(targets, genes_out, exclude=()):
    """Weighted average of centered source changes, on genes_out. Returns (E_c, E_b), 0 where no source."""
    g26 = {g: i for i, g in enumerate(genes26())}
    cols = np.array([g26.get(g, -1) for g in genes_out])
    num_c = np.zeros((len(targets), len(genes_out)))
    num_b, den = np.zeros_like(num_c), np.zeros_like(num_c)
    pos, neg, nsrc = np.zeros_like(num_c), np.zeros_like(num_c), np.zeros_like(num_c)
    for line, w in WEIGHTS.items():
        if line in exclude:
            continue
        c = np.load(OUT / f"{line}.npz")
        rows = {t: i for i, t in enumerate(c["targets"])}
        gw = 1.0
        if EXPR_W is not None:
            z = np.load(ROOT / "data/lines" / f"{line}.npz", allow_pickle=True)
            sc = dict(zip(z["genes"].astype(str), 1e6 * np.asarray(z["ctrl"], float) / np.sum(z["ctrl"])))
            src = np.array([sc.get(g, 0.0) for g in genes_out])
            gw = np.minimum(1.0, 4 * (src + 1) / (EXPR_W[0] + 1)) ** EXPR_W[1]
        for i, t in enumerate(targets):
            if t not in rows:
                continue
            ec = np.where(cols >= 0, c["ec"][rows[t]][cols], np.nan)
            eb = np.where(cols >= 0, c["eb"][rows[t]][cols], np.nan)
            m = ~np.isnan(ec)
            num_c[i] += w * gw * np.where(m, ec, 0)
            pos[i] += m & (ec > 0)
            neg[i] += m & (ec < 0)
            nsrc[i] += m
            num_b[i] += w * gw * np.where(m, eb, 0)
            den[i] += w * gw * m
    with np.errstate(invalid="ignore"):
        ec_f = np.where(den > 0, num_c / den, 0).astype(np.float32)
        eb_f = np.where(den > 0, num_b / den, 0).astype(np.float32)
    fused.agree = np.where(ec_f > 0, pos, neg)  # sources whose sign matches the fused change
    fused.nsrc = nsrc
    return ec_f, eb_f


def promoter_pairs(targets, genes_out):
    """target, neighbor, distance (bp between start sites) for neighbors within 5 kb; written to a CSV."""
    tss = pd.read_csv(ROOT / "data/annot/tss.csv", index_col=0)
    g = tss.reindex(genes_out)
    chrom, pos = g["chr"].to_numpy(), g["tss"].to_numpy()
    rows = []
    for t in targets:
        if t not in tss.index:
            continue
        d = np.abs(pos - tss.loc[t, "tss"])
        near = np.where((chrom == tss.loc[t, "chr"]) & (d <= 5000) & (genes_out != t))[0]
        rows += [(t, genes_out[j], float(d[j])) for j in near]
    path = OUT / f"promoter_pairs_{os.getpid()}.csv"  # per process: parallel runs used to clobber one file
    pd.DataFrame(rows, columns=["target", "neighbor", "distance"]).to_csv(path, index=False)
    return path


def control_stats(X, n_cells, pool, seed):
    """#82's control_template: mean per-cell share m, pooled share q, and n_cells templates that
    each average `pool` random control cells (sorted by depth), with their mean depths."""
    raw = sp.csr_matrix(X, dtype=np.float64)
    lib = np.asarray(raw.sum(1)).ravel()
    m = np.asarray((sp.diags(1 / lib) @ raw).sum(0)).ravel() / len(lib)
    q = np.asarray(raw.sum(0)).ravel()
    q /= q.sum()
    sel = np.random.default_rng(seed).choice(len(lib), n_cells * pool, replace=False)
    sel = sel[np.argsort(lib[sel], kind="stable")]
    tpl = (raw[sel].toarray() / lib[sel, None]).reshape(n_cells, pool, -1).mean(1)
    depths = np.rint(lib[sel].reshape(n_cells, pool).mean(1)).astype(np.int64)
    return m, q, tpl, depths


def profiles(targets, genes_out, m, q, ac, ab, exclude=(), agree=0, thr=0.0):
    """Desired per-cell mean and pooled profile per target (rows sum to 1).
    agree/thr: in the PER-CELL target only, keep a gene's change only if at least `agree` sources measured it
    and agree on its sign, and |change| > thr; the pooled target keeps every change."""
    ec, eb = fused(targets, genes_out, exclude)
    if agree:
        keep = (fused.agree >= agree) & (fused.nsrc >= agree) & (np.abs(ec) > thr)
        gpos = {g: i for i, g in enumerate(genes_out)}
        for i, t in enumerate(targets):  # never mask the knocked-down gene itself
            if t in gpos:
                keep[i, gpos[t]] = True
        ec = np.where(keep, ec, 0).astype(np.float32)
    dc = desired_mean(m, ec, space="log2fc", amplitude=ac, clip=3)
    db = desired_mean(q, eb, space="bulk_delta", amplitude=ab, clip=3)
    pairs = promoter_pairs(targets, genes_out)
    dc, _ = apply_promoter_prior(dc, m, np.asarray(targets), genes_out, pairs, 0.15)
    db, _ = apply_promoter_prior(db, q, np.asarray(targets), genes_out, pairs, 0.15)
    return dc, db


def thin_counts(X, factor, rng):
    """Binomial thinning (Gerard 2020, seqgendiff): real control cells X (cells x genes, counts), each gene's
    counts scaled by `factor` in expectation. factor < 1: keep each count with probability factor;
    factor > 1: add Poisson((factor - 1) * x). Genes with factor 1 stay exactly as the real cells."""
    X = sp.csr_matrix(X, dtype=np.float64)
    f = factor[X.indices]
    x = X.data
    out = np.where(f < 1, rng.binomial(x.astype(np.int64), np.clip(f, 0, 1)),
                   x + rng.poisson(np.maximum(f - 1, 0) * x))
    Y = sp.csr_matrix((out.astype(np.float32), X.indices.copy(), X.indptr.copy()), shape=X.shape)
    Y.eliminate_zeros()
    return Y


def local(line, name, ac=0.6, ab=0.3, pool=4, thin=0, agree=0, thr=0.0):
    """Predict a held-out public line from the others and score it, like src/local_eval.py."""
    import local_eval as le
    out = ROOT / "data/local_eval" / line
    real = ad.read_h5ad(out / "real.h5ad", backed="r")
    genes = np.array(real.var_names, dtype=str)
    counts = real.obs["target"].value_counts()
    targets = [t for t in counts.index if t != le.CTRL]
    ctrl_X = ad.read_h5ad(out / "ctrl_input.h5ad").X
    pool = int(pool)
    n_max = int(counts[targets].max())
    m, q, tpl, depths = control_stats(ctrl_X, n_max, pool, le.SEED)
    dc, db = profiles(targets, genes, m, q, ac, ab, exclude=(line,), agree=int(agree), thr=thr)
    blocks, labels = [sp.csr_matrix(ctrl_X, dtype=np.float32)], [le.CTRL] * ctrl_X.shape[0]
    rng = np.random.default_rng(le.SEED)
    if thin:  # binomial thinning of real control cells by the per-cell target's ratio to the control mean
        Xc = sp.csr_matrix(ctrl_X)
        for i, t in enumerate(targets):
            n = int(counts[t])
            rows = rng.choice(Xc.shape[0], n, replace=False)
            factor = np.divide(dc[i], m, out=np.ones_like(m), where=m > 0)
            blocks.append(thin_counts(Xc[rows], factor, rng))
            labels += [t] * n
    for i, t in enumerate(targets if not thin else []):
        n = int(counts[t])
        idx = np.sort(rng.choice(n_max, n, replace=False)) if n < n_max else np.arange(n_max)
        x = dual_moment_counts(tpl[idx], dc[i], db[i], depths=depths[idx], seed=le.SEED + i)
        blocks.append(sp.csr_matrix(x.astype(np.float32)))
        labels += [t] * n
    pred = ad.AnnData(X=sp.vstack(blocks).tocsr(),
                      obs=pd.DataFrame({"target": labels}, index=np.arange(len(labels)).astype(str)),
                      var=pd.DataFrame(index=genes))
    p = out / f"pred_{name}.h5ad"
    pred.write_h5ad(p)
    print(f"{line}/{name}: ac={ac}, ab={ab}, pool={pool}, thin={thin}, agree={agree}, thr={thr}", flush=True)
    le.cell_eval("run", "-ap", p, "-ar", out / "real.h5ad", "--preset", "vcc2026", "-o", out / f"run_{name}",
                 "--cache-real", out / "real_cache", "--cache-pred", out / f"cache_{name}")
    ref = (["--real-bundle", out / "bundle"] if (out / "bundle").exists() else
           ["--baseline-agg", out / "baseline/baseline_agg.csv",
            "--baseline-meta", out / "baseline/baseline_meta.json"])
    le.cell_eval("score", "--user-agg", out / f"run_{name}" / "agg_results.csv", *ref,
                 "-o", out / f"score_{name}.csv")
    print(pd.read_csv(out / f"score_{name}.csv").to_string())


def write_h5ad64(path, genes, targets, contexts, cells, gen):
    """Slim h5ad like predict_2026.write_h5ad, but CSR with an int64 indptr (nnz can pass 2^31)."""
    import h5py
    n = len(contexts) * len(targets) * cells
    obs = pd.DataFrame({"target_gene": np.tile(np.repeat(targets, cells), len(contexts)),
                        "context": np.repeat(contexts, len(targets) * cells)},
                       index=np.arange(n).astype(str))
    ad.AnnData(X=sp.csr_matrix((n, len(genes)), dtype=np.float32), obs=obs,
               var=pd.DataFrame(index=genes)).write_h5ad(path)
    with h5py.File(path, "r+") as f:
        g = f["X"]
        for k in ("data", "indices", "indptr"):
            del g[k]
        data = g.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
        ind = g.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
        ptr = g.create_dataset("indptr", (n + 1,), dtype="int64")
        ptr[0] = 0
        rows = nnz = 0
        for block in gen:
            b = sp.csr_matrix(block)
            end = nnz + b.nnz
            data.resize((end,))
            ind.resize((end,))
            data[nnz:end] = b.data
            ind[nnz:end] = b.indices
            ptr[rows + 1:rows + b.shape[0] + 1] = b.indptr[1:].astype(np.int64) + nnz
            rows += b.shape[0]
            nnz = end
        assert rows == n
    return n, nnz


def build_2026(name, ac=0.6, ab=0.3, pool=4):
    """The 2026 submission: 300 targets x 400 cells x contexts A, B, C, packaged as .vcc."""
    import predict_2026 as P
    genes = genes26()
    targets = pd.read_csv(ctrl_dir() / "pert_counts.csv")["target_gene"].to_numpy(str)
    pool = int(pool)

    def blocks():
        for k, c in enumerate(P.CONTEXTS):
            a = ad.read_h5ad(P.CTRL_DIR / f"context_{c}.h5ad")
            assert list(a.var_names) == list(genes)
            m, q, tpl, depths = control_stats(a.X, P.CELLS, pool, P.SEED + k)
            dc, db = CTX_PROFILES.get(c, profiles)(list(targets), genes, m, q, ac, ab)
            for i, t in enumerate(targets):
                x = dual_moment_counts(tpl, dc[i], db[i], depths=depths, seed=P.SEED + 1000 * k + i)
                yield sp.csr_matrix(x.astype(np.float32))
                if (i + 1) % 50 == 0:
                    print(f"  context {c}: {i + 1}/{len(targets)} targets", flush=True)

    out = ROOT / "data/submissions"
    h5 = out / f"{name}.h5ad"
    n_obs, nnz = write_h5ad64(h5, genes, targets, P.CONTEXTS, P.CELLS, blocks())
    print(f"wrote {h5} ({n_obs:,} cells, {nnz:,} nonzeros, {nnz / n_obs:,.0f} per cell)", flush=True)
    P.package(h5, out / f"{name}.vcc", n_obs, len(genes), nnz)
    print(f"packaged {out / f'{name}.vcc'}", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build()
    elif sys.argv[1] == "submission":
        opts = dict(a.split("=") for a in sys.argv[3:])
        build_2026(sys.argv[2], **{k: float(v) for k, v in opts.items()})
    else:
        opts = dict(a.split("=") for a in sys.argv[4:])
        local(sys.argv[2], sys.argv[3], **{k: float(v) for k, v in opts.items()})
