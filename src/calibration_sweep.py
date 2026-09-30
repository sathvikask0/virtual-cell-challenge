"""Fast screening, not the official scorer. Shared-gene count and logbulk error.

Each evaluated line is excluded from sources. This is an exploratory tuning
sweep; its best settings require independent confirmation, not a final estimate.
Run: .venv/bin/python src/calibration_sweep.py
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from calibration import mean_prediction
from cross_line import LABELED, fit, load_lines

ROOT = Path(__file__).resolve().parent.parent


def main():
    lines, genes = load_lines(LABELED)
    options = [dict(alpha=a, total_alpha=b, own_alpha=c, mode="transfer")
               for a, b, c in itertools.product([.1, .25, .5, .75, 1.], [0., .25, .5, 1.], [.5, 1.])]
    options += [dict(alpha=.5, total_alpha=.5, own_alpha=1., mode="transfer", expression_gate=g)
                for g in [.1, 1., 5.]]
    options += [dict(mode="control"), dict(mode="target-only")]
    challenge = set(pd.read_csv(ROOT / "data/vcc/controls/pert_counts.csv").target_gene)
    rows = []
    for held, test in lines.items():
        m = fit({n: d for n, d in lines.items() if n != held})
        # Adapt the existing common-gene fit to the production transfer interface.
        src = ({t: np.mean([x for x, _ in v], axis=0) for t, v in m["per_target"].items()},
               {t: np.ones(len(genes)) for t in m["per_target"]},
               {t: [y for _, y in v] for t, v in m["per_target"].items()},
               m["mean_lfc"], m["mean_ltr"], m["own_drop"])
        rng = np.random.default_rng(0)
        seen = [t for t in test["targets"] if t in m["per_target"]]
        groups = {"seen": sorted(rng.choice(seen, min(150, len(seen)), replace=False)),
                  "challenge_overlap": sorted(set(test["targets"]) & challenge)}
        for group, targets in groups.items():
            if not targets:
                continue
            true = test["counts"][[test["tidx"][t] for t in targets]]
            zero_mse = np.square(true - test["ctrl"]).mean()
            for i, opt in enumerate(options):
                pred = np.array([mean_prediction(t, m["gidx"], src, test["ctrl_share"],
                                                 test["ctrl"].sum(), **opt) for t in targets])
                mse = np.square(pred - true).mean()
                # Official v2 comparator, before its sampling correction and gene filter.
                # Keep this explicitly a proxy: scorer applies both to per-cell data.
                logpred = np.log1p(50000 * pred / pred.sum(1, keepdims=True))
                logtrue = np.log1p(50000 * true / true.sum(1, keepdims=True))
                logctrl = np.log1p(50000 * test["ctrl"] / test["ctrl"].sum())
                err = np.square(logpred - logtrue)
                ctrlerr = np.square(logtrue - logctrl)
                for row, target in enumerate(targets):
                    if target in m["gidx"]:
                        err[row, m["gidx"][target]] = 0
                        ctrlerr[row, m["gidx"][target]] = 0
                rows.append(dict(line=held, group=group, setting=i, n=len(targets),
                                 mse=mse, relative_mse=mse / zero_mse,
                                 logbulk_relative_mse=err.sum() / ctrlerr.sum(),
                                 mae=np.abs(pred - true).mean(), **opt))
        print(f"{held}: screened {len(options)} settings", flush=True)
    out = ROOT / "data/calibration"
    out.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "sweep.csv", index=False)
    # Equal weight per cell line; do not let a deep-count line dominate selection.
    print(df[df.group == "seen"].pivot(index="setting", columns="line", values="logbulk_relative_mse")
          .mean(axis=1).sort_values().head(8).to_string())
    (out / "settings.json").write_text(json.dumps(options, indent=2))


if __name__ == "__main__":
    main()
