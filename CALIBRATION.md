# Transfer calibration experiments

Work builds on the existing transfer submission and `cell-eval2` local evaluation.
Original predictors and Claude's denoising work are preserved. New scripts use
`codex_` output names to avoid overwriting existing experiments.

## What is being tested

- Independent shrink factors for expression changes, total RNA, and the targeted gene.
- No-change and target-only comparisons, and soft target-expression gating.
- Control-fitted low-rank Poisson-lognormal noise, compared with independent lognormal
  noise and the existing gamma-Poisson sampler.
- Depth-aware subtraction of Poisson variance when estimating gene overdispersion.

## Scoring correction

The installed `cell-eval2` v0.16.0 scorer's v2 count comparison is
`log1p(50000 * gene_count_sum / total_count_sum)` for MSE and discrimination.
MSE excludes each knockdown's own gene and applies sampling corrections and caps.
PDS excludes all target genes in the panel. DE uses per-cell normalization and
Wilcoxon testing. See the installed package's `metrics/delta.py`,
`metrics/discrimination.py`, and each run's `run_params.yaml`.

Raw-count MSE is useful diagnostically but is not the official expression error.
Changing total RNA alone cannot improve deterministic normalized pseudobulk MSE;
it can still affect sampling noise and detection in generated counts.

## Results so far

The shared-gene sweep evaluates 45 settings on five held-out lines, excluding each
evaluation line from the source profiles. These are exploratory tuning results,
not an independently validated score improvement. Challenge-overlap targets are
reported separately (25 in H1 and 272 in K562).

On HepG2's official 150-target evaluation, changing expression shrinkage from 0.5
to 0.25 while retaining total shrinkage 0.5 reduced overall score from **0.170889
to 0.009691**. PDS stayed close (0.8255 vs 0.8275); the DE metrics deteriorated.
Official normalized MSE also worsened (0.932454 to 0.974504, lower is better).
This candidate is rejected for submission based on this evaluation.

Control-only diagnostic, HepG2, 1,000 held-out control cells, average of three
generation seeds:

| Sampler | Median normalized variance ratio | Mean detected genes | Gene-pair correlation agreement |
|---|---:|---:|---:|
| Real controls | 1 | 4,308 | 1 |
| Original gamma-Poisson | 1.099 | 4,134 | -0.006 |
| Depth-aware gamma-Poisson | 1.034 | 4,242 | -0.002 |
| Independent Poisson-lognormal | 1.097 | 4,160 | 0.000 |
| Factor Poisson-lognormal, strength 0.5 | 1.098 | 4,182 | 0.285 |

Pair agreement is Pearson correlation of generated versus real gene-pair
correlations on 200 genes chosen using input controls only. Variance ratios use
genes with held-out normalized mean above 1 per 10,000 counts. This diagnoses
control generation, not perturbation accuracy. Neither noise variant is a
submission winner until the full official score confirms it.

The H1 control diagnostic confirms both changes generalize beyond HepG2:
depth-aware noise improves the variance ratio from 1.041 to 0.994 and detected
genes from 8,681 to 8,751 (real: 8,742). Factor noise improves gene-pair agreement
from -0.005 to 0.256, but does not fix marginal variance by itself.

Full HepG2 factor-noise and depth-aware-noise evaluations are in progress.
No new leaderboard upload yet. Nine behavior and simulation tests pass, including
recovery of dispersion under variable depth and the correlated sampler's expected
means, variances, and covariance.

## Reproduce

Run from the repository root:

```sh
.venv/bin/python -m unittest discover -s src -p 'test_*.py'
.venv/bin/python src/calibration_sweep.py
.venv/bin/python src/control_diagnostics.py hepg2
.venv/bin/python src/control_diagnostics.py h1
.venv/bin/python src/calibration_eval.py hepg2 a025 --alpha .25 --total-alpha .5
.venv/bin/python src/calibration_eval.py hepg2 factor_a05 --alpha .5 --sampler factor
.venv/bin/python src/calibration_eval.py hepg2 depth_a05 --alpha .5 --depth-aware
```

Prediction files are never silently overwritten. To rerun scoring on an existing
prediction, append `--score-only`. Add `--bulk-only` for official MSE and PDS
without expensive DE or an overall score. The latter changes only the requested
metric subset; it is not used to claim an overall improvement.

Per-run settings are saved as JSON beside predictions. Detailed results and
diagnostics live under `data/local_eval/` and `data/calibration/` (gitignored).
Source profiles use only other public cell lines; factors and variance fits use
only the selected line's input controls. Strong candidates need confirmation
across multiple seeds and held-out cell lines before leaderboard upload.
