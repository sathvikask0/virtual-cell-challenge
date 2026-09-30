"""Official-scoring companion to calibration_sweep; separate from shared scripts.

Example: python src/calibration_eval.py hepg2 a025 --alpha .25 --total-alpha .5
Reuse an existing prediction with --score-only. Names never overwrite old files.
"""
import argparse
import gc
import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from calibration import calibrated_change
from correlated_noise import fit_factors, sample_lognormal
from local_eval import ROOT, CTRL, cell_eval
from predict_2026 import EPS, SOURCES, context_stats, sample_cells, source_changes
from sampling_stats import depth_aware_stats
from prediction_io import write_prediction_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("line", choices=["h1", "hepg2", "jurkat"])
    parser.add_argument("name")
    parser.add_argument("--alpha", type=float, default=.25)
    parser.add_argument("--total-alpha", type=float, default=.5)
    parser.add_argument("--own-alpha", type=float, default=1.)
    parser.add_argument("--expression-gate", type=float, default=0.)
    parser.add_argument("--mode", choices=["transfer", "control", "target-only"], default="transfer")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--sampler", choices=["gamma", "lognormal", "factor"], default="gamma")
    parser.add_argument("--factor-strength", type=float, default=.5)
    parser.add_argument("--depth-aware", action="store_true")
    parser.add_argument("--effect-space", choices=["original", "logbulk"], default="original")
    parser.add_argument("--k562-source", choices=["original", "core", "batch"], default="original")
    parser.add_argument("--score-only", action="store_true")
    parser.add_argument("--bulk-only", action="store_true",
                        help="Fast official MSE/PDS only; no overall score")
    args = parser.parse_args()
    if args.effect_space == "logbulk" and (args.own_alpha != 1. or args.expression_gate != 0.):
        parser.error("logbulk currently transfers the source's own-gene effect; own-alpha and expression-gate are supported only in original space")
    out = ROOT / "data/local_eval" / args.line
    path = out / f"pred_codex_{args.name}.h5ad"
    if not args.score_only:
        if path.exists():
            raise FileExistsError(f"Use a new experiment name or --score-only: {path}")
        real = ad.read_h5ad(out / "real.h5ad", backed="r")
        genes = np.array(real.var_names, dtype=str)
        counts = real.obs["target"].value_counts()
        real.file.close()
        targets = [t for t in counts.index if t != CTRL]
        sources = [s for s in SOURCES if s != args.line]
        if args.k562_source != "original":
            sources = [f"k562_{args.k562_source}" if s == "k562" else s for s in sources]
        src = source_changes(genes, sources)
        logbulk_effects = {}
        if args.effect_space == "logbulk":
            # Transfer differences in the expression comparator used by the scorer.
            # Include only requested targets and retain source measurement masks.
            wanted = set(targets)
            for source in sources:
                with np.load(ROOT / "data/lines" / f"{source}.npz") as data:
                    positions = {g: i for i, g in enumerate(genes)}
                    mask = np.array([g in positions for g in data["genes"]])
                    cols = np.array([positions[g] for g in data["genes"][mask]])
                    control = data["ctrl"][mask]
                    control_log = np.log1p(50000 * control / control.sum())
                    for row, target in enumerate(data["targets"]):
                        if target not in wanted:
                            continue
                        counts_row = data["counts"][row, mask]
                        delta = np.log1p(50000 * counts_row / counts_row.sum()) - control_log
                        if target not in logbulk_effects:
                            logbulk_effects[target] = (np.zeros(len(genes)), np.zeros(len(genes)))
                        summed, measured = logbulk_effects[target]
                        summed[cols] += delta
                        measured[cols] += 1
        gidx = {g: i for i, g in enumerate(genes)}
        ctrl = ad.read_h5ad(out / "ctrl_input.h5ad").X
        share, totals, phi = (depth_aware_stats if args.depth_aware else context_stats)(ctrl)
        factors = fit_factors(ctrl, phi) if args.sampler == "factor" else None
        rng = np.random.default_rng(args.seed)
        labels = [CTRL] * ctrl.shape[0]
        for t in targets:
            labels.extend([t] * int(counts[t]))
        opts = {k: getattr(args, k) for k in
                ["alpha", "total_alpha", "own_alpha", "expression_gate", "mode"]}
        def blocks():
            yield sp.csr_matrix(ctrl, dtype=np.float32)
            for t in targets:
                lfc, ltr = calibrated_change(t, gidx, src, share, **opts)
                if args.effect_space == "logbulk" and t in logbulk_effects and args.mode == "transfer":
                    summed, measured = logbulk_effects[t]
                    delta = summed / np.maximum(measured, 1)
                    s = np.expm1(np.maximum(np.log1p(50000 * share) + args.alpha * delta, 0))
                else:
                    s = np.maximum((share + EPS) * np.exp(lfc) - EPS, 0)
                if args.sampler == "gamma":
                    yield sample_cells(s / s.sum(), totals, phi, ltr, rng, n=int(counts[t]))
                else:
                    yield sample_lognormal(s / s.sum(), totals, phi, ltr, rng, n=int(counts[t]),
                                          factors=factors, strength=args.factor_strength)
        write_prediction_blocks(path, genes, labels, blocks())
        (out / f"config_codex_{args.name}.json").write_text(json.dumps(vars(args), indent=2))
        print(f"Wrote {path.name}: {opts}", flush=True)
        # Do not retain multi-GB source profiles while the scorer subprocess runs.
        del src, blocks, ctrl, factors, logbulk_effects
        gc.collect()
    run = out / f"run_codex_{args.name}{'_bulk' if args.bulk_only else ''}"
    extra = (["--set", "metrics=[pds_cosine,expr_mse_unbiased_capped_norm]"]
             if args.bulk_only else [])
    cell_eval("run", "-ap", path, "-ar", out / "real.h5ad", "--preset", "vcc2026", "-o", run, *extra)
    if args.bulk_only:
        print(pd.read_csv(run / "agg_results.csv").to_string(index=False))
        return
    ref = (["--real-bundle", out / "bundle"] if (out / "bundle").exists() else
           ["--baseline-agg", out / "baseline/baseline_agg.csv",
            "--baseline-meta", out / "baseline/baseline_meta.json"])
    score = out / f"score_codex_{args.name}.csv"
    cell_eval("score", "--user-agg", run / "agg_results.csv", *ref, "-o", score)
    print(pd.read_csv(score).to_string(index=False))


if __name__ == "__main__":
    main()
