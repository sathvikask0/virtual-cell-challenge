"""Learn how a knockdown's change carries over from one cell line to another (gradient boosting).

One row = (destination line D, target t, gene g). Features come from the other big lines
(K562, HCT116, HEK293T minus D): their change for (t, g), its z-score, gene g's normal expression
in the sources and in D, target t's expression in D, own-gene and neighbor flags. Label = D's
change for (t, g), centered (each line's average change over all its targets removed), because
the part shared by all targets doesn't help tell targets apart.

Training targets exclude the 2026 targets and the H1 targets, so the H1 check is clean.

Commands:
  python src/gbm.py build     cache per-line changes for the targets used (data/gbm/{line}.npz)
  python src/gbm.py train     fit, then the fast H1 check vs plain K562
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/gbm"
EPS = 1e-5
BIG = ["k562", "hct116", "hek293t"]   # lines that cover the 2026 targets (sources)
DEST = BIG + ["rpe1", "hepg2", "jurkat"]
N_TRAIN = 1500                        # training targets
GENES_PER_TARGET = 400                # random genes per (D, t), plus every gene a source flags
SEED = 0


def genes26():
    return pd.read_csv(ROOT / "data/vcc/controls/gene_names.csv")["gene_name"].to_numpy(str)


def build():
    OUT.mkdir(exist_ok=True)
    genes = genes26()
    gidx = {g: i for i, g in enumerate(genes)}
    t26 = set(pd.read_csv(ROOT / "data/vcc/controls/pert_counts.csv")["target_gene"])
    th1 = set(np.load(ROOT / "data/lines/h1.npz")["targets"])
    big_targets = [set(np.load(ROOT / "data/lines" / f"{l}.npz")["targets"]) for l in BIG]
    # training targets: knocked down in at least 2 big lines, not 2026, not H1
    cand = sorted(t for t in set.union(*big_targets)
                  if sum(t in s for s in big_targets) >= 2 and t not in t26 and t not in th1)
    rng = np.random.default_rng(SEED)
    train = set(rng.choice(cand, N_TRAIN, replace=False))
    keep = train | t26 | th1
    for line in DEST + ["h1"]:
        d = np.load(ROOT / "data/lines" / f"{line}.npz")
        cols = np.array([gidx.get(g, -1) for g in d["genes"]])
        ok = cols >= 0
        counts, ctrl = d["counts"][:, ok], d["ctrl"][ok]
        share = counts / counts.sum(1, keepdims=True)
        q = ctrl / ctrl.sum()
        lfc = np.log((share + EPS) / (q + EPS))
        mean = lfc.mean(0)
        nz = np.load(ROOT / "data/lines" / f"{line}_noise.npz")
        nidx = {g: i for i, g in enumerate(nz["genes"])}
        var = nz["var"][[nidx[g] for g in d["genes"][ok]]]
        rows = np.array([t in keep for t in d["targets"]])
        r = np.log(counts[rows].sum(1) / ctrl.sum())
        diff = counts[rows] - ctrl[None] * np.exp(r)[:, None]
        z = diff / np.sqrt(np.maximum(var, 1e-6)[None] * (1 / d["n"][rows, None] + 1 / d["ctrl_n"]))

        def full(a):  # onto the 2026 gene list, NaN where the line doesn't measure the gene
            out = np.full(a.shape[:-1] + (len(genes),), np.nan, np.float32)
            out[..., cols[ok]] = a
            return out

        np.savez(OUT / f"{line}.npz", targets=d["targets"][rows],
                 lfc=full(lfc[rows] - mean), z=full(z), mean=full(mean),
                 logq=full(np.log(q + EPS)), train=np.array(sorted(train & set(d["targets"]))))
        print(f"{line}: {rows.sum()} targets cached", flush=True)
        del d, counts, share, lfc


def load():
    return {l: dict(np.load(OUT / f"{l}.npz")) for l in DEST + ["h1"]}


def neighbor_dist(genes):
    tss = pd.read_csv(ROOT / "data/annot/tss.csv", index_col=0).reindex(genes)
    return tss["chr"].to_numpy(), tss["tss"].to_numpy()


def features(L, dest, t, gsel, srcs, chrom, pos, genes):
    """Feature matrix for target t, genes gsel, destination dest, from source lines srcs.
    Aggregated over sources, so the same model serves any set of sources."""
    lf, zz, lq, mn = [], [], [], []
    for s in srcs:
        r = L[s]["row"].get(t)
        lf.append(L[s]["lfc"][r, gsel] if r is not None else np.full(len(gsel), np.nan))
        zz.append(L[s]["z"][r, gsel] if r is not None else np.full(len(gsel), np.nan))
        lq.append(L[s]["logq"][gsel])
        mn.append(L[s]["mean"][gsel])
    lf, zz, lq, mn = map(np.array, (lf, zz, lq, mn))
    f = {}
    with np.errstate(all="ignore"):
        f["lfc_avg"] = np.nanmean(lf, 0)
        f["lfc_sd"] = np.nanstd(lf, 0)
        best = np.nanargmax(np.where(np.isnan(zz), -1, np.abs(zz)), 0)
        f["lfc_best"] = lf[best, np.arange(len(gsel))]
        f["z_best"] = zz[best, np.arange(len(gsel))]
        f["z_avg"] = np.nanmean(zz, 0)
        f["n_src"] = np.sum(~np.isnan(lf), 0).astype(np.float32)
        f["logq_src"] = np.nanmean(lq, 0)
        f["mean_src"] = np.nanmean(mn, 0)
    f["logq_dest"] = L[dest]["logq"][gsel]
    ti = L["_gidx"].get(t)
    f["target_logq_dest"] = np.full(len(gsel), L[dest]["logq"][ti] if ti is not None else np.nan)
    with np.errstate(all="ignore"):
        f["target_logq_src"] = np.full(len(gsel), np.nanmean([L[s]["logq"][ti] for s in srcs])
                                       if ti is not None else np.nan)
    f["is_own"] = (genes[gsel] == t).astype(np.float32)
    if ti is not None and isinstance(chrom[ti], str):
        d = np.where(chrom[gsel] == chrom[ti], np.abs(pos[gsel] - pos[ti]), np.inf)
    else:
        d = np.full(len(gsel), np.inf)
    f["log_dist"] = np.log10(np.minimum(d, 1e8) + 1).astype(np.float32)
    return pd.DataFrame(f)


def train():
    genes = genes26()
    L = load()
    for l in L:
        L[l]["row"] = {t: i for i, t in enumerate(L[l]["targets"])}
    L["_gidx"] = {g: i for i, g in enumerate(genes)}
    chrom, pos = neighbor_dist(genes)
    rng = np.random.default_rng(SEED)
    Xs, ys, ws, grp = [], [], [], []
    train_t = set(L["k562"]["train"]) | set(L["hct116"]["train"]) | set(L["hek293t"]["train"])
    for dest in DEST:
        srcs = [s for s in BIG if s != dest]
        measured = np.where(~np.isnan(L[dest]["logq"]))[0]
        for t in L[dest]["targets"]:
            if t not in train_t:
                continue
            flagged = np.zeros(len(genes), bool)
            for s in srcs:
                r = L[s]["row"].get(t)
                if r is not None:
                    flagged |= np.abs(np.nan_to_num(L[s]["z"][r])) > 3
            flagged &= ~np.isnan(L[dest]["logq"])
            rand = rng.choice(measured, GENES_PER_TARGET, replace=False)
            gsel = np.unique(np.concatenate([np.where(flagged)[0], rand]))
            w = np.where(flagged[gsel], 1.0, len(measured) / GENES_PER_TARGET)
            y = L[dest]["lfc"][L[dest]["row"][t], gsel]
            m = ~np.isnan(y)
            Xs.append(features(L, dest, t, gsel[m], srcs, chrom, pos, genes))
            ys.append(np.clip(y[m], -3, 3))
            ws.append(w[m])
            grp.append(np.full(m.sum(), DEST.index(dest)))
        print(f"rows so far {sum(map(len, ys)):,} (after {dest})", flush=True)
    X = pd.concat(Xs, ignore_index=True).astype(np.float32)
    y, w = np.concatenate(ys), np.concatenate(ws)
    del Xs, ys, ws
    model = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, max_leaf_nodes=63,
                                          min_samples_leaf=200, early_stopping=True,
                                          validation_fraction=0.1, random_state=SEED)
    model.fit(X, y, sample_weight=w)
    print(f"fit: {model.n_iter_} trees, {len(y):,} rows", flush=True)
    import joblib
    joblib.dump((model, list(X.columns)), OUT / "model.joblib")
    del X, y, w
    check_h1(model, L, genes, chrom, pos)


def predict(model, L, dest, t, gsel, chrom, pos, genes):
    """Average the model over pairs of big source lines (training saw at most 2 sources)."""
    pairs = [[a, b] for i, a in enumerate(BIG) for b in BIG[i + 1:] if dest not in (a, b)]
    out = [model.predict(features(L, dest, t, gsel, p, chrom, pos, genes).astype(np.float32))
           for p in pairs]
    return np.mean(out, 0)


def changes_for(targets, genes_out, ctrl_share, scale_k562=False):
    """Model-predicted change for each target on genes_out, for a new context whose control
    share (on genes_out) is ctrl_share. Targets no big line covers are left out."""
    import joblib
    model, _ = joblib.load(OUT / "model.joblib")
    genes = genes26()
    L = load()
    for l in L:
        L[l]["row"] = {t: i for i, t in enumerate(L[l]["targets"])}
    L["_gidx"] = {g: i for i, g in enumerate(genes)}
    oidx = {g: i for i, g in enumerate(genes_out)}
    logq = np.full(len(genes), np.nan, np.float32)
    have = np.array([g in oidx for g in genes])
    logq[have] = np.log(ctrl_share[[oidx[g] for g in genes[have]]] + EPS)
    L["_ctx"] = {"logq": logq}
    chrom, pos = neighbor_dist(genes)
    gsel = np.where(have)[0]
    cols = np.array([oidx[g] for g in genes[gsel]])
    out = {}
    for t in targets:
        if not any(t in L[s]["row"] for s in BIG):
            continue
        g = predict(model, L, "_ctx", t, gsel, chrom, pos, genes)
        if scale_k562 and t in L["k562"]["row"]:
            k = np.nan_to_num(L["k562"]["lfc"][L["k562"]["row"][t], gsel])
            g = g * np.linalg.norm(k) / max(np.linalg.norm(g), 1e-9)
        v = np.zeros(len(genes_out))
        v[cols] = g
        out[t] = v
    return out


def check_h1(model, L, genes, chrom, pos):
    """Fast pds-like check on H1: share of other targets whose true profile is closer (lower is
    better), vs the plain centered K562 change. H1 targets were never used in training."""
    h = np.load(ROOT / "data/lines/h1.npz")
    hidx = {g: i for i, g in enumerate(h["genes"])}
    gsel = np.array([i for i, g in enumerate(genes) if g in hidx])
    hc = np.array([hidx[g] for g in genes[gsel]])
    q = h["ctrl"][hc] / h["ctrl"][hc].sum()
    hs = h["counts"][:, hc] / h["counts"][:, hc].sum(1, keepdims=True)
    rh = {t: i for i, t in enumerate(h["targets"])}
    te = [t for t in h["targets"] if t in L["k562"]["row"]]
    f = lambda p: np.log1p(1e4 * p)
    gpos = {g: i for i, g in enumerate(genes[gsel])}
    preds = {"K562 (current)": [], "GBM": [], "GBM, K562 size": []}
    truth = []
    for t in te:
        k = np.nan_to_num(L["k562"]["lfc"][L["k562"]["row"][t], gsel])
        g = predict(model, L, "h1", t, gsel, chrom, pos, genes)
        scale = np.linalg.norm(k) / max(np.linalg.norm(g), 1e-9)
        for name, l in zip(preds, [k, g, g * scale]):
            s = np.clip((q + EPS) * np.exp(l) - EPS, 0, None)
            p = f(s / s.sum()) - f(q)
            if t in gpos:
                p[gpos[t]] = 0
            preds[name].append(p)
        y = f(hs[rh[t]]) - f(q)
        if t in gpos:
            y[gpos[t]] = 0
        truth.append(y)
    T = np.array(truth)
    n = lambda A: A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    rng = np.random.default_rng(0)
    base = None
    for name, P in preds.items():
        S = n(np.array(P)) @ n(T).T
        v = np.array([(S[j] > S[j, j]).mean() for j in range(len(te))])
        msg = f"{name:16s} share closer {v.mean():.3f}  cos {np.diag(S).mean():.3f}"
        if base is None:
            base = v
        else:
            d = v - base
            bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
            msg += f"  vs K562 {d.mean():+.3f} [{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]"
        print(msg, flush=True)


if __name__ == "__main__":
    {"build": build, "train": train}[sys.argv[1]]()
