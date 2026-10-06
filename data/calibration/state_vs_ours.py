"""Head-to-head on held-out Jurkat: Arc State (st-se-replogle-full/jurkat_0.99 published predictions) vs our recipe,
on the targets both cover, scored against our real Jurkat local eval. Run from repo root:
  .venv/bin/python data/calibration/state_vs_ours.py [our_run]
Metrics per method: target-specific correlation (panel-mean removed), retrieval among the common targets,
sign accuracy on real-significant genes ranked by |predicted LFC| (top 50), split by DepMap essentiality."""
import glob
import sys

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

RUN = sys.argv[1] if len(sys.argv) > 1 else "rmpc5"
real = pl.read_parquet(glob.glob("data/local_eval/jurkat/real_cache/de_wilcoxon_table-*.parquet")[0])
ours = pl.read_parquet(glob.glob(f"data/local_eval/jurkat/cache_{RUN}/de_wilcoxon_table-*.parquet")[0])

a = ad.read_h5ad("data/state/jurkat_pred.h5ad", backed="r")
pcol = "gene" if "gene" in a.obs else a.obs.columns[0]
lab = a.obs[pcol].astype(str).to_numpy()
import pickle
genes_s = np.array(pickle.load(open("data/state/var_dims.pkl", "rb"))["gene_names"], str)
T = sorted((set(real["target"].unique().to_list()) & set(lab)) - {"non-targeting"})
print(f"State obs col {pcol!r}, {a.n_obs} cells, {len(genes_s)} genes; common targets {len(T)}")


def mean_expr(mask):
    idx = np.flatnonzero(mask)
    acc = np.zeros(a.n_vars)
    for s0 in range(0, len(idx), 5000):
        X = a.X[idx[s0:s0 + 5000]]
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        acc += np.expm1(X).sum(0)  # log1p space -> linear
    return acc / max(len(idx), 1)


print("non-targeting cells in State pred:", (lab == "non-targeting").sum())
ctrl = mean_expr(lab == "non-targeting")
S = {}
for t in T:
    m = mean_expr(lab == t)
    S[t] = np.log2((m + 1e-3) / (ctrl + 1e-3))
S = pd.DataFrame(S, index=genes_s).T


def table(df, col="log2_fold_change"):
    w = df.filter(pl.col("target").is_in(T)).select("target", "feature", col).to_pandas()
    return w.pivot(index="target", columns="feature", values=col)


R = table(real)
O = table(ours)
sig = table(real, "p_adj")
G = sorted(set(R.columns) & set(O.columns) & set(S.columns))
G = [g for g in G if g not in T]
R, O, S, sig = R.loc[T, G].fillna(0), O.loc[T, G].fillna(0), S.loc[T, G].fillna(0), sig.loc[T, G].fillna(1)
print(f"genes compared {len(G)}")

GE = pd.read_csv("data/depmap/CRISPRGeneEffect.csv", index_col=0)
GE.columns = [c.split(" (")[0] for c in GE.columns]
ess = GE.mean().reindex(T)


def evaluate(P):
    Pc, Rc = P - P.mean(0), R - R.mean(0)
    corr = np.array([np.corrcoef(Pc.iloc[i], Rc.iloc[i])[0, 1] for i in range(len(T))])
    Pn = Pc.values / (np.linalg.norm(Pc.values, axis=1, keepdims=True) + 1e-12)
    Rn = Rc.values / (np.linalg.norm(Rc.values, axis=1, keepdims=True) + 1e-12)
    Sim = Pn @ Rn.T
    ret = np.array([1 - (Sim[i] > Sim[i, i]).sum() / len(T) for i in range(len(T))])
    acc = []
    for i in range(len(T)):
        m = sig.values[i] < .05
        if m.sum() < 5:
            acc.append(np.nan)
            continue
        idx = np.flatnonzero(m)
        top = idx[np.argsort(-np.abs(P.values[i, idx]))[:50]]
        acc.append(np.mean(np.sign(P.values[i, top]) == np.sign(R.values[i, top])))
    return pd.DataFrame({"corr": corr, "retrieval": ret, "sign_top50": acc}, index=T)


res = {"ours": evaluate(O), "State": evaluate(S)}
for lab_, sel in [("all", ess.notna() | ess.isna()), ("non-essential (> -.3)", ess > -.3), ("essential (<= -.3)", ess <= -.3)]:
    print(f"== {lab_}: {int(sel.sum())} targets")
    for k, v in res.items():
        d = v[sel.to_numpy()]
        print(f"   {k:6s} corr {d['corr'].mean():.3f}  retrieval {d.retrieval.mean():.3f}  sign_top50 {d.sign_top50.mean():.3f}")
