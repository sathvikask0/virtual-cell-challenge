"""Package a calibrated candidate using the existing streaming .vcc writer.

Example: python src/calibration_submission.py codex_depth_a1 --depth-aware
No upload is performed. Use --pilot-targets 2 for a small archive/schema check.
"""
import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd

from calibration import calibrated_change
from correlated_noise import fit_factors, sample_lognormal
from predict_2026 import (ROOT, CTRL_DIR, CONTEXTS, EPS, context_stats, source_changes,
                          sample_cells, write_h5ad, package)
from sampling_stats import depth_aware_stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("--alpha", type=float, default=1.)
    parser.add_argument("--total-alpha", type=float, default=1.)
    parser.add_argument("--own-alpha", type=float, default=1.)
    parser.add_argument("--sampler", choices=["gamma", "lognormal", "factor"], default="gamma")
    parser.add_argument("--factor-strength", type=float, default=.5)
    parser.add_argument("--depth-aware", action="store_true")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--pilot-targets", type=int, default=0)
    args = parser.parse_args()
    out = ROOT / "data/submissions"
    out.mkdir(exist_ok=True)
    h5, archive = out / f"{args.name}.h5ad", out / f"{args.name}.vcc"
    if h5.exists() or archive.exists():
        raise FileExistsError("Choose a new candidate name; existing artifacts are preserved")
    genes = pd.read_csv(CTRL_DIR / "gene_names.csv").gene_name.to_numpy(str)
    targets = pd.read_csv(CTRL_DIR / "pert_counts.csv").target_gene.to_numpy(str)
    if args.pilot_targets < 0 or args.pilot_targets > len(targets):
        parser.error("pilot-targets must be between 0 and the full target count")
    if args.pilot_targets:
        targets = targets[:args.pilot_targets]
    src = source_changes(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(args.seed)

    def blocks():
        for context in CONTEXTS:
            controls = ad.read_h5ad(CTRL_DIR / f"context_{context}.h5ad")
            if list(controls.var_names) != list(genes):
                raise ValueError(f"Gene order differs in context {context}")
            stats = depth_aware_stats if args.depth_aware else context_stats
            share, totals, phi = stats(controls.X)
            factors = fit_factors(controls.X, phi) if args.sampler == "factor" else None
            for i, target in enumerate(targets):
                lfc, ltr = calibrated_change(target, gidx, src, share, alpha=args.alpha,
                                            total_alpha=args.total_alpha, own_alpha=args.own_alpha)
                s = np.maximum((share + EPS) * np.exp(lfc) - EPS, 0)
                s /= s.sum()
                if args.sampler == "gamma":
                    yield sample_cells(s, totals, phi, ltr, rng)
                else:
                    yield sample_lognormal(s, totals, phi, ltr, rng, factors=factors,
                                           strength=args.factor_strength)
                if (i + 1) % 50 == 0 or i + 1 == len(targets):
                    print(f"Context {context}: {i + 1}/{len(targets)} targets", flush=True)

    config = out / f"{args.name}.json"
    config.write_text(json.dumps(vars(args), indent=2))
    n_obs, nnz = write_h5ad(h5, genes, targets, blocks())
    package(h5, archive, n_obs, len(genes), nnz)
    config.write_text(json.dumps(dict(vars(args), n_obs=n_obs, n_genes=len(genes), nnz=nnz,
                                      complete=True, pilot=bool(args.pilot_targets)), indent=2))
    print(f"Packaged {archive}: {n_obs:,} cells, {len(genes):,} genes", flush=True)


if __name__ == "__main__":
    main()
