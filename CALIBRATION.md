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

Full official HepG2 results at expression/total shrinkage 0.5:

| Sampler | Overall | PDS scaled | Raw normalized MSE (lower better) |
|---|---:|---:|---:|
| Original gamma-Poisson | 0.170889 | 0.825507 | 0.932454 |
| Factor Poisson-lognormal | 0.182908 | 0.853001 | 0.916680 |
| Depth-aware gamma-Poisson | 0.184269 | 0.902273 | 0.874819 |

The improvements are about the same size as the previously observed seed
variation (0.015); they are promising, not confirmed submission improvements.
An independent lognormal bulk-only run gives MSE 0.921736 and raw PDS 0.784027,
versus factor lognormal 0.916680 and 0.781029. Therefore the correlated variant's
small bulk gains cannot all be attributed to its correlations.

Concurrent experiments from Claude found alpha=1 gives a much larger official
HepG2 gain: 0.373523. The next comparisons apply depth-aware noise at alpha=1,
total-alpha=1, and repeat the original alpha=1 model with a second seed.

Those comparisons have now finished: original alpha=1 scores 0.373523 (seed 0)
and 0.389126 (seed 1), whereas depth-aware alpha=1 scores 0.360574 (seed 0).
Depth-aware noise is **not selected** for the next submission: matching control
moments better did not improve the stronger model's full score. The original
no-shrink model is the current baseline. A small gain at alpha=0.5 is insufficient
evidence to replace it.

User reports Claude submitted the next candidate and requests no more submissions
today (2026-09-30). Remaining work is local only. New H1 bulk-only comparisons
reuse the existing predictions without duplicating the expensive DE computation.

H1 official bulk-only results: alpha=0.5 gives raw normalized MSE 1.378686 and
raw PDS 0.826890; alpha=1 gives MSE 2.705667 and PDS 0.830112. The full H1
baseline-relative scores later completed at 0.176660 and 0.178456 respectively.
These are not on HepG2's replicate scale, so do not compare their numeric levels
across lines. Alpha=2 drops the H1 baseline-relative score to 0.069161.

## Further local experiments

`learned_transfer.py` fits a 32-component effect basis and a ridge correction
using ordered source/destination pairs for the same knockdown. Each outer test
line is excluded from both PCA and regression fitting. The context variant has
133 input features; the ablation retains 34 (source effect and target expression).
Both variants hurt discrimination on H1, HepG2, and Jurkat. At correction strength
0.5, the context model's proxy PDS changes from direct-copy 0.895 to 0.857 on H1,
0.764 to 0.687 on HepG2, and 0.818 to 0.790 on Jurkat. Mean error improves slightly
on the cancer lines, but these models are not selected for submission. All these
are shared-gene, uncorrected mean-expression proxies, not overall scores.

An official HepG2 evaluation of `--effect-space logbulk` is in progress. It copies
differences in `log1p(50000 * share)` directly rather than EPS-regularized share
ratios. In this variant the target gene's own effect is also copied from sources;
independent own-gene scaling and expression gating are currently unsupported.

The authors' normalized K562 pseudobulk file (file 35773217) was downloaded and
its MD5 verified as `a3dfaa94ea8724217f5ecb1e14a5f0c8`. It has the same 8,248 raw
gene columns as the raw bulk file (8,246 unique symbols); the 66 GB single-cell
file is not needed to access this processed alternative.

`build_k562_sources.py` writes separate sources, preserving `k562.npz`:

- `k562_core.npz`: original target averages, using the authors' 514 vetted control
  guides, excluding 71 other control guides.
- `k562_batch.npz`: authors' gemgroup-normalized effects mapped to a count proxy
  with scales fitted to those vetted controls. This conversion is experimental,
  not observed counts. The 6,881 nonfinite normalized entries are treated as
  unavailable effect estimates, not extreme gene changes.

Source provenance and processing diagnostics are saved in
`data/calibration/k562_sources.json`. H1 and HepG2 batch-source bulk-only
evaluations are in progress. The local prediction writer now streams cell blocks
to HDF5, avoiding retention and duplication of entire prediction matrices.
Thirteen behavior and numerical tests pass.

The new `calibration_submission.py` reuses the original streaming packager and
supports calibrated shrinkage, all samplers, and depth-aware noise. A two-target,
three-context pilot passed `vcc prep --dry-run` (2,400 integer-count cells, 18,533
genes, correct context coverage and gene order). This pilot disables verification
of the complete official target list; a full candidate must enable that check.

No new leaderboard upload by Codex yet. Nine behavior and simulation tests pass, including
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
.venv/bin/python src/calibration_eval.py hepg2 depth_a1 --alpha 1 --total-alpha 1 --depth-aware
.venv/bin/python src/calibration_eval.py hepg2 gamma_a1_s1 --alpha 1 --total-alpha 1 --seed 1
.venv/bin/python src/calibration_submission.py codex_depth_pilot --depth-aware --pilot-targets 2
.venv/bin/vcc prep data/submissions/codex_depth_pilot.h5ad --dry-run --no-verify-targets --genes data/vcc/controls/gene_names.csv
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
