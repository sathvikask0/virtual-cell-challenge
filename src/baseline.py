"""Score simple baselines on held-out knocked-down genes, in counts space.

Every prediction is built from two parts:
  share  what fraction of the cell's RNA each measured gene is  (18080 numbers, sum to 1)
  total  average total counts per cell                          (1 number)
  predicted counts per cell = share * total

The true answer for gene G is the average raw counts per cell of G's cells.

Baselines (fit on train genes only, scored on val and test genes):
  control:  predict the control cells (knockdown changes nothing)
  mean:     control share + average train share change,
            control total * average train total ratio

Metrics, averaged over held-out genes:
  pearson_delta   correlation of predicted vs true change from control, in counts
                  (1 = right shape of response, 0 = no relation)
  mae             mean absolute error of predicted vs true counts per cell
  total_err       |predicted total / true total - 1|, e.g. 0.05 = total off by 5%
  discrimination  rank of the true row among all held-out rows by L1 distance,
                  scaled to 0..1 (0 = perfect, 0.5 = random guessing)
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CTRL = "non-targeting"


def load(split):
    """Return targets, share, mean_total, mean_counts for knockdowns, plus the same for control."""
    d = np.load(ROOT / "data/pseudobulk" / f"{split}.npz")
    targets = d["targets"]
    share = d["sum_counts"] / d["sum_total"][:, None]
    total = d["sum_total"] / d["n"]
    counts = d["sum_counts"] / d["n"][:, None]
    c, k = targets == CTRL, targets != CTRL
    ctrl = dict(share=share[c][0], total=total[c][0], counts=counts[c][0])
    return targets[k], share[k], total[k], counts[k], ctrl


def score(pred_counts, true_counts, ctrl_counts):
    p = pred_counts - ctrl_counts
    p[np.abs(p) < 1e-9] = 0  # rounding noise is not a predicted change
    t = true_counts - ctrl_counts
    p, t = p - p.mean(1, keepdims=True), t - t.mean(1, keepdims=True)
    denom = np.linalg.norm(p, axis=1) * np.linalg.norm(t, axis=1)
    pearson = np.where(denom > 0, (p * t).sum(1) / np.maximum(denom, 1e-12), 0.0)

    mae = np.abs(pred_counts - true_counts).mean(1)
    total_err = np.abs(pred_counts.sum(1) / true_counts.sum(1) - 1)

    # dist[i, j] = L1 distance between prediction i and true row j (in chunks to save memory)
    dist = np.concatenate([np.abs(pred_counts[i:i + 16, None, :] - true_counts[None]).sum(2)
                           for i in range(0, len(pred_counts), 16)])
    own = np.diag(dist)
    # ties count as half, so identical predictions for everyone score 0.5
    rank = (dist < own[:, None]).sum(1) + 0.5 * ((dist == own[:, None]).sum(1) - 1)
    disc = rank / (len(own) - 1)
    return dict(pearson_delta=pearson.mean(), mae=mae.mean(),
                total_err=total_err.mean(), discrimination=disc.mean())


def main():
    _, tr_share, tr_total, _, tr_ctrl = load("train")
    share_change = (tr_share - tr_ctrl["share"]).mean(0)
    total_ratio = (tr_total / tr_ctrl["total"]).mean()

    def control(n, ctrl):
        return np.repeat(ctrl["counts"][None], n, axis=0)

    def mean(n, ctrl):
        share = ctrl["share"] + share_change
        return np.repeat((share * ctrl["total"] * total_ratio)[None], n, axis=0)

    print(f"train: average total ratio (knockdown / control) = {total_ratio:.3f}\n")
    print(f"{'split':6} {'baseline':9} {'pearson_delta':>14} {'mae':>8} {'total_err':>10} {'discrimination':>15}")
    for split in ["val", "test"]:
        targets, _, _, true_counts, ctrl = load(split)
        for name, make in {"control": control, "mean": mean}.items():
            s = score(make(len(targets), ctrl), true_counts, ctrl["counts"])
            print(f"{split:6} {name:9} {s['pearson_delta']:14.3f} {s['mae']:8.3f} "
                  f"{s['total_err']:10.3f} {s['discrimination']:15.3f}")


if __name__ == "__main__":
    main()
