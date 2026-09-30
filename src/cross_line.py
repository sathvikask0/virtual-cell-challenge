"""Practice the 2026 task: predict knockdowns in a cell line the model has never seen.

Leave one cell line out: fit on the other lines, predict the held-out line using only its
control cells, score against its true answers. Only genes measured in every line are used.

Changes are learned as ratios of share (fraction of the cell's RNA), because lines differ
~5x in counts per cell:
  lfc[g]      = log((knockdown share[g] + EPS) / (control share[g] + EPS))
  total ratio = log(knockdown total / control total)
Prediction = held-out control share * exp(lfc), renormalized, * control total * exp(total ratio).

Predictors:
  control     no change
  mean        average lfc over all knockdowns (each training line weighted equally)
  mean+drop   mean, plus the switched-off gene itself drops by the typical amount
  transfer    if this gene was switched off in a training line, reuse its average lfc there;
              otherwise mean+drop

Held-out targets are split into those switched off in some training line ("seen") and
those that never were ("unseen"). Up to MAX_EVAL random targets of each are scored.

Usage: python src/cross_line.py [line ...]   (default: every line with answers in data/lines)
"""
import sys
from pathlib import Path

import numpy as np

from baseline import score

ROOT = Path(__file__).resolve().parent.parent
LABELED = ["h1", "k562", "rpe1", "hepg2", "jurkat"]
EPS = 1e-5
MAX_EVAL = 300


def load_lines(names):
    lines = {n: dict(np.load(ROOT / "data/lines" / f"{n}.npz")) for n in names}
    genes = set.intersection(*(set(d["genes"]) for d in lines.values()))
    for n in ["ctx_A", "ctx_B", "ctx_C"]:  # keep only genes the 2026 data also measures
        p = ROOT / "data/lines" / f"{n}.npz"
        if p.exists():
            genes &= set(np.load(p)["genes"])
    genes = np.array(sorted(genes))
    for d in lines.values():
        idx = {g: i for i, g in enumerate(d["genes"])}
        cols = [idx[g] for g in genes]
        d["counts"], d["ctrl"] = d["counts"][:, cols], d["ctrl"][cols]
        d["genes"] = genes
        share = d["counts"] / d["counts"].sum(1, keepdims=True)
        d["ctrl_share"] = d["ctrl"] / d["ctrl"].sum()
        d["lfc"] = np.log((share + EPS) / (d["ctrl_share"] + EPS))
        d["ltr"] = np.log(d["counts"].sum(1) / d["ctrl"].sum())
        d["tidx"] = {t: i for i, t in enumerate(d["targets"])}
    return lines, genes


def fit(train):
    """Everything the predictors need, from the training lines only."""
    gidx = {g: i for i, g in enumerate(next(iter(train.values()))["genes"])}
    mean_lfc = np.mean([d["lfc"].mean(0) for d in train.values()], axis=0)
    mean_ltr = np.mean([d["ltr"].mean() for d in train.values()])
    own = [d["lfc"][i, gidx[t]] for d in train.values()
           for t, i in d["tidx"].items() if t in gidx]
    per_target = {}
    for d in train.values():
        for t, i in d["tidx"].items():
            per_target.setdefault(t, []).append((d["lfc"][i], d["ltr"][i]))
    return dict(gidx=gidx, mean_lfc=mean_lfc, mean_ltr=mean_ltr,
                own_drop=np.median(own), per_target=per_target)


def predict(kind, m, targets, ctrl_share, ctrl_total):
    rows = []
    for t in targets:
        lfc, ltr = m["mean_lfc"].copy(), m["mean_ltr"]
        if kind == "control":
            lfc, ltr = np.zeros_like(lfc), 0.0
        elif kind == "transfer" and t in m["per_target"]:
            lfc = np.mean([x for x, _ in m["per_target"][t]], axis=0)
            ltr = np.mean([y for _, y in m["per_target"][t]])
        elif kind in ("mean+drop", "transfer") and t in m["gidx"]:
            lfc[m["gidx"][t]] = m["own_drop"]
        share = np.clip((ctrl_share + EPS) * np.exp(lfc) - EPS, 0, None)
        rows.append(share / share.sum() * ctrl_total * np.exp(ltr))
    return np.array(rows)


def main():
    names = sys.argv[1:] or [n for n in LABELED if (ROOT / "data/lines" / f"{n}.npz").exists()]
    lines, genes = load_lines(names)
    rng = np.random.default_rng(0)
    print(f"lines: {', '.join(names)}   shared genes: {len(genes):,}\n")
    print(f"{'held out':9} {'targets':12} {'predictor':10} {'pearson_delta':>14} {'mae':>7} "
          f"{'total_err':>10} {'discrim':>8}")
    for held in names:
        test = lines[held]
        m = fit({n: d for n, d in lines.items() if n != held})
        seen = [t for t in test["targets"] if t in m["per_target"]]
        unseen = [t for t in test["targets"] if t not in m["per_target"]]
        for label, group in [("seen", seen), ("unseen", unseen)]:
            if len(group) < 10:
                continue
            group = list(rng.choice(group, min(MAX_EVAL, len(group)), replace=False))
            true = test["counts"][[test["tidx"][t] for t in group]]
            for kind in ["control", "mean", "mean+drop", "transfer"]:
                if label == "unseen" and kind == "transfer":
                    continue  # same as mean+drop
                pred = predict(kind, m, group, test["ctrl_share"], test["ctrl"].sum())
                s = score(pred, true, test["ctrl"])
                print(f"{held:9} {label + f' ({len(group)})':12} {kind:10} "
                      f"{s['pearson_delta']:14.3f} {s['mae']:7.3f} {s['total_err']:10.3f} "
                      f"{s['discrimination']:8.3f}", flush=True)
        print()


if __name__ == "__main__":
    main()
