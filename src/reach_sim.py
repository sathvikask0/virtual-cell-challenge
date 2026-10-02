"""Offline reach (cell_eval2 de_direction_reach, raw, universe='adjudicated') from cached DE tables, and what-if
re-rankings of the predicted |log2FC|.
Per target: the REAL significant genes (p_adj < .05, own gene excluded), ranked by (our p_adj < .05 first, then our
|LFC| descending); purity = running share of matching signs; k* = deepest prefix with purity >= .9; reach = k*/N_conf.
Usage (from repo root): python src/reach_sim.py LINE RUN
"""
import glob
import sys

import anndata as ad
import numpy as np
import polars as pl


def tables(line, run):
    p = pl.read_parquet(glob.glob(f"data/local_eval/{line}/cache_{run}/de_wilcoxon_table-*.parquet")[0])
    r = pl.read_parquet(glob.glob(f"data/local_eval/{line}/real_cache/de_wilcoxon_table-*.parquet")[0])
    return p, r


def reach(p, r, key=None):
    """key: optional polars frame (target, feature, key) replacing |LFC| as the ranking score."""
    real = r.filter((pl.col("p_adj") < .05) & (pl.col("target") != pl.col("feature")) & (pl.col("log2_fold_change") != 0))
    real = real.select("target", "feature", pl.col("log2_fold_change").alias("lr"))
    pr = p.select("target", "feature", pl.col("log2_fold_change").alias("lp"), pl.col("p_adj").alias("pp"))
    d = real.join(pr, on=["target", "feature"], how="left")
    if key is not None:
        d = d.join(key, on=["target", "feature"], how="left")
    else:
        d = d.with_columns(pl.col("lp").abs().alias("key"))
    d = d.with_columns(((pl.col("pp") < .05).fill_null(False)).alias("sig"),
                       (np.sign(pl.col("lp")) == np.sign(pl.col("lr"))).fill_null(False).alias("ok"),
                       pl.col("key").fill_null(-1.0))
    out = []
    for (t,), g in d.group_by(["target"]):
        g = g.sort(["sig", "key"], descending=[True, True])
        ok = g["ok"].to_numpy().astype(float)
        pur = np.cumsum(ok) / np.arange(1, len(ok) + 1)
        good = np.flatnonzero(pur >= .9)
        out.append((good[-1] + 1 if len(good) else 0) / len(ok))
    return float(np.mean(out)), len(out)


if __name__ == "__main__":
    line, run = sys.argv[1], sys.argv[2]
    p, r = tables(line, run)
    print(f"{line}/{run}: offline raw reach {reach(p, r)[0]:.4f}")
