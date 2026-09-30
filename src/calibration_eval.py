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
    parser.add_argument("--score-only", action="store_true")
    parser.add_argument("--bulk-only", action="store_true",
                        help="Fast official MSE/PDS only; no overall score")
    args = parser.parse_args()
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
        src = source_changes(genes, [s for s in SOURCES if s != args.line])
        gidx = {g: i for i, g in enumerate(genes)}
        ctrl = ad.read_h5ad(out / "ctrl_input.h5ad").X
        share, totals, phi = (depth_aware_stats if args.depth_aware else context_stats)(ctrl)
        factors = fit_factors(ctrl, phi) if args.sampler == "factor" else None
        rng = np.random.default_rng(args.seed)
        blocks, labels = [sp.csr_matrix(ctrl, dtype=np.float32)], [CTRL] * ctrl.shape[0]
        opts = {k: getattr(args, k) for k in
                ["alpha", "total_alpha", "own_alpha", "expression_gate", "mode"]}
        for t in targets:
            lfc, ltr = calibrated_change(t, gidx, src, share, **opts)
            s = np.maximum((share + EPS) * np.exp(lfc) - EPS, 0)
            if args.sampler == "gamma":
                block = sample_cells(s / s.sum(), totals, phi, ltr, rng, n=int(counts[t]))
            else:
                block = sample_lognormal(s / s.sum(), totals, phi, ltr, rng, n=int(counts[t]),
                                        factors=factors, strength=args.factor_strength)
            blocks.append(block)
            labels.extend([t] * int(counts[t]))
        pred = ad.AnnData(X=sp.vstack(blocks).tocsr(),
                          obs=pd.DataFrame({"target": labels}, index=np.arange(len(labels)).astype(str)),
                          var=pd.DataFrame(index=genes))
        pred.write_h5ad(path)
        (out / f"config_codex_{args.name}.json").write_text(json.dumps(vars(args), indent=2))
        print(f"Wrote {path.name}: {opts}", flush=True)
        # Do not retain multi-GB source profiles while the scorer subprocess runs.
        del src, pred, blocks, ctrl, factors
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
