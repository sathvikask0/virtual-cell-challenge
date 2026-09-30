"""Compare generated unchanged controls to independent real controls.

This diagnoses the existing sampler without perturbation-prediction error.
Run: python src/control_diagnostics.py hepg2
"""
import argparse
import json

import anndata as ad
import numpy as np
import scipy.sparse as sp

from local_eval import ROOT, CTRL
from predict_2026 import context_stats, sample_cells
from correlated_noise import fit_factors, sample_lognormal
from sampling_stats import depth_aware_stats


def summaries(X, cols):
    X = sp.csr_matrix(X, dtype=np.float64)
    totals = np.asarray(X.sum(1)).ravel()
    norm = sp.diags(10000 / np.maximum(totals, 1)) @ X
    mean = np.asarray(norm.mean(0)).ravel()
    variance = np.asarray(norm.multiply(norm).mean(0)).ravel() - mean ** 2
    selected = norm[:, cols].toarray()
    correlation = np.corrcoef(selected, rowvar=False)
    offdiag = correlation[np.triu_indices(len(cols), 1)]
    return dict(n_cells=X.shape[0], median_total=float(np.median(totals)),
                mean_detected=float(np.diff(X.indptr).mean()),
                median_abs_gene_correlation=float(np.nanmedian(np.abs(offdiag)))), mean, variance, correlation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("line", choices=["h1", "hepg2", "jurkat"])
    args = parser.parse_args()
    out = ROOT / "data/local_eval" / args.line
    ctrl = ad.read_h5ad(out / "ctrl_input.h5ad").X
    share, totals, phi = context_stats(ctrl)
    _, _, adjusted_phi = depth_aware_stats(ctrl)
    real = ad.read_h5ad(out / "real.h5ad", backed="r")
    indices = np.flatnonzero(real.obs.target.astype(str).to_numpy() == CTRL)
    rng = np.random.default_rng(123)
    indices = np.sort(rng.choice(indices, min(1000, len(indices)), replace=False))
    reference = sp.csr_matrix(real.X[indices])
    real.file.close()
    # Select features from model input only; avoid constant/very sparse columns.
    candidates = np.flatnonzero(share * 10000 > 1)
    cols = np.sort(rng.choice(candidates, min(200, len(candidates)), replace=False))
    ref, refmean, refvar, refcorr = summaries(reference, cols)
    factors = fit_factors(ctrl, phi)
    result = {"line": args.line, "real": ref}
    for method in ["gamma_poisson", "depth_aware_gamma", "independent_lognormal", "factor_lognormal"]:
        result[method] = []
        for seed in [0, 1, 2]:
            generator_rng = np.random.default_rng(seed)
            if method in ("gamma_poisson", "depth_aware_gamma"):
                p = adjusted_phi if method == "depth_aware_gamma" else phi
                generated = sample_cells(share, totals, p, 0., generator_rng, n=len(indices))
            else:
                generated = sample_lognormal(share, totals, phi, 0., generator_rng, n=len(indices),
                                            factors=factors if method == "factor_lognormal" else None,
                                            strength=.5)
            summary, mean, var, corr = summaries(generated, cols)
            mask = refmean > 1
            upper = np.triu_indices(len(cols), 1)
            summary.update(seed=seed, median_variance_ratio=float(np.median(var[mask] / refvar[mask])),
                           gene_pair_correlation_agreement=float(np.corrcoef(refcorr[upper], corr[upper])[0, 1]),
                           normalized_mean_mae=float(np.abs(mean - refmean).mean()))
            result[method].append(summary)
    destination = ROOT / "data/calibration"
    destination.mkdir(exist_ok=True)
    (destination / f"controls_{args.line}.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
