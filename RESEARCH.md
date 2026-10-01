# Research notes (2026-10-01)

What the literature and the scorer's code say about where our remaining score is.
Current leaderboard: 0.159 (pds 0.64, mse 0.07, nmae 0.16, fid −0.02, reach 0.11, jac 0.00).

## Papers and write-ups

- **Molina & Zhang 2026, "Perturbation response decomposition…"** (bioRxiv 10.64898/2026.07.24.740459), four CRISPRi
  screens (K562, RPE1, HepG2, Jurkat). They split the reproducible response into four parts:

  | part | share | can a zero-shot model predict it? |
  |---|---|---|
  | shared by all knockdowns within a cell line | 28% | no; the scorer's baseline already has the true one |
  | specific to the target, same in every cell line | 29% | **yes, this is our ceiling** |
  | target × cell line interaction | 24% | no, not without data from that line |
  | noise | 19% | no |

  Simple ridge/MLP models beat complex ones. Pooling 3 source lines beats 1. The new line's control expression
  doesn't help (matches our tests). Essential-gene knockdowns are about 2x more consistent across lines.
  → More, cleaner sources for the shared target-specific part is the main lever. That's Codex's work (CD4 etc.).
- **Ahlmann-Eltze et al. 2025, Nature Methods:** deep models don't beat simple linear baselines. That matches our
  GBM / ridge / PCA results.
- **"Locked evaluation… sampling-depth entanglement"** (arXiv 2608.00152): the VCC-style measure is strongly tied to
  cell count and depth, and in-distribution success didn't transfer. A warning: significance-based metrics react to
  how noisy and deep the predicted cells are, not only to the means.
- **VCC 2025 wrap-up (Arc):** almost all models were worse than baseline on MAE; top teams used classical statistics
  plus some deep learning; "global linear scaling" was used to optimize PDS.
- **Public #82 (0.155):** we adopted their method → 0.159.
- **#1 Illumina AI (0.42):** no public method. Most likely large proprietary Perturb-seq data. Not something we can copy.

## How the significant-gene metrics work (cell_eval2/metrics/direction.py, de.py)

- Significance comes from a Wilcoxon test of predicted (or real) cells vs real controls, Benjamini–Hochberg,
  p_adj < 0.05, with the target's own gene excluded.
- **fid (direction_fidelity_yield_raw)** = k / max(n_pred, n_conf):
  - n_pred = genes we call significant;
  - k = those whose sign matches the real log2FC;
  - n_conf = genes significant in the real data.
  → Calling far more genes than the truth has is penalized directly.
- **reach** ranks genes by our significance and finds the deepest prefix that is still ≥ 90% sign-correct, as a share of n_conf.
- **jac** = overlap of the predicted and real significant sets.
- All three are scaled so that the true-mean-response baseline = 0.

## Measured on H1 (true DE tables saved by Codex)

- The true number of significant genes per target is very uneven: median 268, mean 1,238, 10th pct 4, 90th pct 4,294.
- The old method called a median of 782 genes for every target. For weak targets that's far too many, so fid collapses.
- The shared response is **not** what fills the true lists: only 2% of true significant genes are significant in
  ≥ 25% of targets. Adding a shared response to the per-cell target has little room.

## Next experiment

**Match how many genes we make significant to each target's expected strength.** Strong knockdowns (big, consistent
source changes) should get many significant genes; weak ones few. First a diagnostic: rerun the current best
(ac 1.0 / ab 0.5) on H1 with the DE tables saved (`--cache-pred`), and compare our per-target count with the true count
and with the strength of the source signal. If we over-call weak targets, scale each target's per-cell change
(or add noise to its cells) by its predicted strength.
