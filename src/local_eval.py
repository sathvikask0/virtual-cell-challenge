"""Score predictions locally with Arc's official scorer (cell-eval2, preset vcc2026).

A held-out cell line plays the role of a 2026 context:
  - pick N_TARGETS random targets, up to CELLS_PER_TARGET real cells each
  - split its control cells in two: one half is the model's input (like the 2026 download),
    the other half goes with the answers (the scorer uses held-out controls, as Arc does)
  - the model predicts from the other public lines only (the held-out line is excluded)

Commands:
  python src/local_eval.py setup LINE [--bundle]
                                            write data/local_eval/LINE/{real.h5ad, ctrl_input.h5ad},
                                            build the official baseline (0 end of the scale) and,
                                            with --bundle, the replicate (1 end; >2 h on a laptop)
  python src/local_eval.py predict LINE NAME [alpha=0.5] [k=K] [kt=KT] [gamma=1] [phi_scale=1]
                                            (k, kt, gamma: options of predict_2026.py;
                                             phi_scale: multiply the sampler's overdispersion)
                                            write a transfer prediction and score it

LINE: h1 (VCC 2025, per-cell) or hepg2 / jurkat (Nadig 2025, per-cell).
"""
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from predict_2026 import (SOURCES, context_stats, recenter, sample_cells, source_changes,
                          target_change)

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / ".venv/bin"
CTRL = "non-targeting"
N_TARGETS = 150
CELLS_PER_TARGET = 200
N_CTRL = 5000
SEED = int(os.environ.get("SEED", 0))  # sampling seed for predictions (setup always uses 0)
FILES = {
    "h1": (ROOT / "data/vcc/train/adata_Training.h5ad", "target_gene"),
    "hepg2": (ROOT / "data/nadig/hepg2_raw_singlecell.h5ad", "gene"),
    "jurkat": (ROOT / "data/nadig/jurkat_raw_singlecell.h5ad", "gene"),
}


def cell_eval(*args):
    cmd = [str(BIN / "cell-eval2"), *map(str, args)]
    print("  $", " ".join(cmd[1:3]), "...", flush=True)
    subprocess.run(cmd, check=True)


def setup(line):
    out = ROOT / "data/local_eval" / line
    out.mkdir(parents=True, exist_ok=True)
    path, col = FILES[line]
    a = ad.read_h5ad(path, backed="r")
    labels = a.obs[col].astype(str).to_numpy()
    genes = np.array(a.var["gene_name"] if "gene_name" in a.var else a.var_names, dtype=str)
    _, first = np.unique(genes, return_index=True)
    keep_genes = np.sort(first)  # drop duplicate gene symbols

    rng = np.random.default_rng(0)
    targets = rng.choice(sorted(set(labels) - {CTRL}), N_TARGETS, replace=False)
    rows = []
    for t in targets:
        idx = np.where(labels == t)[0]
        rows.append(np.sort(rng.choice(idx, min(CELLS_PER_TARGET, len(idx)), replace=False)))
    ctrl = rng.permutation(np.where(labels == CTRL)[0])
    n = min(N_CTRL, len(ctrl) // 2)
    ctrl_real, ctrl_input = np.sort(ctrl[:n]), np.sort(ctrl[n:2 * n])

    def take(idx):
        X = sp.csr_matrix(a.X[np.sort(idx)], dtype=np.float32)
        return X[:, keep_genes]

    real_idx = np.concatenate(rows + [ctrl_real])
    order = np.argsort(real_idx)
    X = take(real_idx)  # rows in sorted index order
    lab = labels[real_idx][order]
    real = ad.AnnData(X=X, obs=pd.DataFrame({"target": lab}, index=np.arange(len(lab)).astype(str)),
                      var=pd.DataFrame(index=genes[keep_genes]))
    real.write_h5ad(out / "real.h5ad")
    ad.AnnData(X=take(ctrl_input), var=pd.DataFrame(index=genes[keep_genes])).write_h5ad(
        out / "ctrl_input.h5ad")
    print(f"{line}: {len(targets)} targets, {X.shape[0]:,} real cells ({n:,} controls), "
          f"{n:,} input controls, {len(keep_genes):,} genes")

    cell_eval("baseline", "-ar", out / "real.h5ad", "--preset", "vcc2026",
              "-o", out / "baseline", "--save-pred", out / "baseline_pred.h5ad")
    if "--bundle" not in sys.argv:  # the replicate end of the scale; slow (5 extra full runs)
        return
    cell_eval("prep-real-bundle", "--real", out / "real.h5ad",
              "--baseline", out / "baseline_pred.h5ad", "-o", out / "bundle",
              "--preset", "vcc2026", "--force")


def predict(line, name, alpha=0.5, k=None, kt=None, gamma=1.0, phi_scale=1.0):
    out = ROOT / "data/local_eval" / line
    real = ad.read_h5ad(out / "real.h5ad", backed="r")
    genes = np.array(real.var_names, dtype=str)
    counts = real.obs["target"].value_counts()
    targets = [t for t in counts.index if t != CTRL]
    src = source_changes(genes, [s for s in SOURCES if s != line], k, kt)
    gidx = {g: i for i, g in enumerate(genes)}
    ctrl_X = ad.read_h5ad(out / "ctrl_input.h5ad").X
    share, totals, phi = context_stats(ctrl_X)
    phi = phi * phi_scale  # <1: less cell-to-cell spread than the controls
    eps = 1e-5
    rng = np.random.default_rng(SEED)
    # the scorer needs control rows in the prediction file too, but scores against the
    # real held-out controls (control_source: real), so these input controls are not used
    blocks, labels = [sp.csr_matrix(ctrl_X, dtype=np.float32)], [CTRL] * ctrl_X.shape[0]
    changes = {t: target_change(t, gidx, src, alpha) for t in targets}
    changes = recenter(changes, gidx, gamma)
    for t in targets:
        lfc, r = changes[t]
        s = np.clip((share + eps) * np.exp(lfc) - eps, 0, None)
        blocks.append(sample_cells(s / s.sum(), totals, phi, r, rng, n=int(counts[t])))
        labels += [t] * int(counts[t])
    seen = sum(t in src[0] for t in targets)
    pred = ad.AnnData(X=sp.vstack(blocks).tocsr(),
                      obs=pd.DataFrame({"target": labels}, index=np.arange(len(labels)).astype(str)),
                      var=pd.DataFrame(index=genes))
    p = out / f"pred_{name}.h5ad"
    pred.write_h5ad(p)
    print(f"{line}/{name}: alpha={alpha}, k={k}, kt={kt}, {seen}/{len(targets)} targets seen in other lines")
    cell_eval("run", "-ap", p, "-ar", out / "real.h5ad", "--preset", "vcc2026",
              "-o", out / f"run_{name}")
    # with the full bundle: official scale (baseline 0, replicate 1); without: baseline only
    ref = (["--real-bundle", out / "bundle"] if (out / "bundle").exists() else
           ["--baseline-agg", out / "baseline/baseline_agg.csv",
            "--baseline-meta", out / "baseline/baseline_meta.json"])
    cell_eval("score", "--user-agg", out / f"run_{name}" / "agg_results.csv", *ref,
              "-o", out / f"score_{name}.csv")
    print(pd.read_csv(out / f"score_{name}.csv").to_string())


if __name__ == "__main__":
    cmd, line = sys.argv[1], sys.argv[2]
    if cmd == "setup":
        setup(line)
    else:
        opts = dict(a.split("=") for a in sys.argv[4:])
        predict(line, sys.argv[3], **{key: float(v) for key, v in opts.items()})
