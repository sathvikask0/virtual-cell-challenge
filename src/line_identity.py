"""Layer 0: place each context in DepMap (CCLE) expression space from its control pseudobulk.

Per context: log1p(CPM) of summed control counts, on genes shared with DepMap 24Q4 log2(TPM+1).
Both sides are centred per gene by the DepMap mean, restricted to the most variable DepMap genes,
then Pearson-correlated against every DepMap line. Known lines (Jurkat, HepG2, K562, HCT116, HEK293T)
are a sanity check: they should map to themselves.

Usage (from src/): python line_identity.py [n_genes] [name ...]   (names default to A/B/C + known lines)
"""
import sys

import numpy as np
import pandas as pd

import atlas_shift as A

D = A.ROOT / "data/depmap"
NAMES = ["ctx_A", "ctx_B", "ctx_C", "h1", "jurkat", "hepg2", "k562", "hct116", "hek293t", "rpe1"]


def depmap():
    X = pd.read_csv(D / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
    X = X.loc[:, [c for c in X.columns if "(" in c]] if any("(" in c for c in X.columns) else X
    X.columns = [c.split(" (")[0] for c in X.columns]
    X = X.loc[:, ~X.columns.duplicated()]
    M = pd.read_csv(D / "Model.csv", index_col=0)
    return X, M


def profile(name):
    z = np.load(A.ROOT / "data/lines" / f"{name}.npz", allow_pickle=True)
    c = np.asarray(z["ctrl"], float)
    return pd.Series(np.log2(1e6 * c / c.sum() + 1), index=z["genes"].astype(str))


def main(n_genes=2000, names=NAMES):
    X, M = depmap()
    label = M.reindex(X.index)
    for name in names:
        try:
            p = profile(name)
        except FileNotFoundError:
            continue
        p = p[~p.index.duplicated()]
        g = X.columns.intersection(p.index)
        Xg = X[g]
        top = Xg.var().sort_values(ascending=False).index[:n_genes]
        mu = Xg[top].mean()
        a = (p[top] - mu).to_numpy()
        B = (Xg[top] - mu).to_numpy()
        r = (B @ a) / (np.linalg.norm(B, axis=1) * np.linalg.norm(a) + 1e-12)
        o = np.argsort(-r)[:5]
        hits = "; ".join(f"{label['CellLineName'].iloc[i]} ({label['OncotreeLineage'].iloc[i]}, "
                         f"{label['Sex'].iloc[i]}) {r[i]:.2f}" for i in o)
        lin = pd.Series(r, index=label["OncotreeLineage"].fillna("?")).groupby(level=0).max().sort_values(ascending=False)
        print(f"{name:8s} genes {len(g)}  top: {hits}")
        print(f"{'':8s} best per lineage: " + ", ".join(f"{k} {v:.2f}" for k, v in lin.head(4).items()), flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2000, sys.argv[2:] or NAMES)
