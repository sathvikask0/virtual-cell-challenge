# Experiments

Task: predict how every gene's counts change when one gene is switched off (CRISPRi),
for genes the model has never seen switched off.

Data: Virtual Cell Challenge 2025 (H1 stem cells), `data/vcc/`. Scripts: `src/`.

| split | genes switched off | cells | use |
|---|---|---|---|
| train | 150 | 221,273 | fit models, cross-validation |
| val | 50 | 98,927 | compare ideas |
| test | 100 | 170,846 | final check only |

The same 38,176 control cells (nothing switched off) appear in all three files.
Gene sets do not overlap between splits.

## How we score

Everything is in **counts space**: predicted counts per cell = share × total.
- **share**: fraction of the cell's RNA each of the 18,080 genes makes up
- **total**: average total counts per cell

The answer key for gene G is the average raw counts per cell of G's cells (= share × total exactly).

| metric | meaning | better |
|---|---|---|
| pearson_delta | correlation of predicted vs true change from control | higher (1 = perfect) |
| mae | average error per gene, in counts per cell | lower |
| total_err | how far off the total counts per cell is (0.04 = 4%) | lower |
| discrimination | can you tell which gene was switched off from the prediction? | lower (0 = perfect, 0.5 = random) |

---

## Done

### Exp 1: mean baseline, log space (superseded by Exp 2)
- Each cell scaled to 10k total, log1p, averaged per switched-off gene, minus control average.
- Prediction: average change across the 150 train genes, same for every gene.
- Result: pearson_delta 0.170 (val) / 0.174 (test), discrimination 0.500.
- Sanity check: the switched-off gene itself drops in all 50 val genes (median −0.68 log units ≈ half).
- Dropped because the 2026 scorer reportedly works in raw counts, not log.

### Finding: total counts per cell
- Control cells vary ~2.4× in total counts (10th–90th pct: 35k–84k), mostly measurement noise.
- Lab batch alone shifts the median total 43k–65k.
- 80% of switched-off genes change the total by less than −8% / +3%.
- A few lower it 20–30%: PRDM14 (0.70), KDM1A (0.74), METTL14 (0.78), SUPT4H1 (0.79), METTL3 (0.79).
- Decision: model share and total separately, multiply to get counts.

### Exp 2: mean baseline, counts space (share × total) — the bar to beat
- `src/pseudobulk.py`, `src/baseline.py`
- share = control share + average train share change; total = control total × average train total ratio (0.979).

| split | predictor | pearson_delta | mae | total_err | discrimination |
|---|---|---|---|---|---|
| val | control (no change) | 0.000 | 0.194 | 0.036 | 0.500 |
| val | **mean** | **0.145** | 0.207 | 0.039 | 0.500 |
| test | control (no change) | 0.000 | 0.194 | 0.040 | 0.500 |
| test | **mean** | **0.159** | 0.205 | 0.042 | 0.500 |

### Exp 3: co-expression embedding + kNN / ridge
- `src/embed_coexpr.py`: gene embedding from control cells only. Scale to 10k, log1p, z-score each gene,
  truncated SVD (128 dims). Genes that rise and fall together get nearby vectors.
- `src/embed_models.py`: input = embedding of the switched-off gene; output = share change + log total ratio.
  - knn: average response of the k most similar train genes
  - ridge: linear regression
- Settings (dims 16/64/128, k 1–20, alpha 0.1–100) picked by 5-fold CV over train genes, then scored on val. Test untouched.
- Note: the "mean" row here is 0.141 on val, not 0.145, because it averages the log total ratio and clips/renormalizes share.

5-fold CV on 150 train genes (best rows):

| model | setting | pearson_delta | mae | total_err | discrimination |
|---|---|---|---|---|---|
| mean | — | 0.217 | 0.209 | 0.041 | 0.498 |
| knn | k=20, 64 dims | 0.221 | 0.219 | 0.043 | 0.501 |
| ridge | alpha=10, 128 dims (64 ties) | 0.232 | 0.210 | 0.041 | 0.492 |
| ridge | alpha=0.1, 128 dims (lowest discrim) | 0.176 | 0.290 | 0.048 | 0.439 |

Val (50 held-out genes), best CV setting per model:

| model | setting | pearson_delta | mae | total_err | discrimination |
|---|---|---|---|---|---|
| mean | — | 0.141 | 0.209 | 0.039 | 0.500 |
| knn | k=20, 64 dims | 0.152 | 0.217 | 0.040 | 0.502 |
| **ridge** | alpha=10, 128 dims | **0.168** | **0.206** | **0.036** | **0.490** |

- Result: ridge beats the mean on every metric, but by a small margin. kNN barely helps.
- Discrimination is still close to random (0.49). Settings that tell genes apart better (low alpha, small k)
  predict each gene's change worse, so predictions stay close to the average.

### Finding: the 2026 task is different from our H1 setup
- 2026 has no training set. Validation = 3 hidden cell lines (A, B, C) in `data/vcc/controls/`:
  18,400 control cells each (raw integer counts, 18,533 genes) + the same 300 genes to switch off in all three.
  Final test = 3 other cell lines (D, E, F), released Oct 22.
- Only 25 of the 300 targets were switched off in H1 (13 train, 4 val, 8 test). All 300 are measured genes in H1.
- 18,077 of the 18,533 genes overlap with H1's 18,080.
- Median total counts per control cell is ~20k in all three lines (10th–90th pct ~9k–34k), vs ~50k in H1.
- So our H1 val (new genes, same cell line) does not test the real task (mostly new genes, new cell line).

### Exp 4: practice on unseen cell lines (leave one cell line out)
- Goal: copy the 2026 setup. Fit on the other cell lines, predict the held-out one from its control cells only.
- Data (public, free), `src/lines.py` puts each in one format in `data/lines/`:

| line | source | targets | control counts/cell |
|---|---|---|---|
| h1 | VCC 2025 | 300 | 57.5k |
| k562 | Replogle 2022 genome-wide (averaged file) | 9,866 | 11.7k |
| rpe1 | Replogle 2022 essential (averaged file) | 2,393 | 13.0k |
| hepg2 | Nadig 2025 essential (per-cell file) | 2,393 | 19.2k |
| jurkat | Nadig 2025 essential (per-cell file) | 2,393 | 12.5k |

- 272 of the 300 2026 targets were switched off in K562. None were in the essential screens.
- `src/cross_line.py`: scored on the 6,123 genes measured in every line (Replogle files keep only ~8k genes).
  Changes are learned as ratios of share, since counts per cell differ ~5x between lines.
- Predictors: control (no change), mean (average change), mean+drop (switched-off gene itself drops),
  transfer (reuse the same gene's average change from the training lines; else mean+drop).
- "seen" = the gene was switched off in some training line. Up to 300 random targets scored per group.

| held out | targets | predictor | pearson_delta | mae | discrimination |
|---|---|---|---|---|---|
| h1 | seen (280) | control | 0.00 | **0.406** | 0.500 |
| h1 | seen (280) | mean+drop | 0.12 | 0.542 | 0.495 |
| h1 | seen (280) | transfer | **0.22** | 0.719 | **0.361** |
| k562 | seen (300) | control | 0.00 | **0.130** | 0.500 |
| k562 | seen (300) | mean+drop | 0.26 | 0.135 | 0.498 |
| k562 | seen (300) | transfer | **0.27** | 0.171 | **0.339** |
| rpe1 | seen (300) | control | 0.00 | 0.285 | 0.500 |
| rpe1 | seen (300) | mean+drop | 0.03 | 0.282 | 0.499 |
| rpe1 | seen (300) | transfer | **0.24** | **0.277** | **0.409** |
| hepg2 | seen (300) | control | 0.00 | 0.330 | 0.500 |
| hepg2 | seen (300) | mean+drop | 0.17 | **0.318** | 0.499 |
| hepg2 | seen (300) | transfer | **0.32** | 0.332 | **0.343** |
| jurkat | seen (300) | control | 0.00 | 0.236 | 0.500 |
| jurkat | seen (300) | mean+drop | 0.11 | 0.236 | 0.499 |
| jurkat | seen (300) | transfer | **0.39** | 0.239 | **0.339** |
| k562 | unseen (300) | control | 0.00 | **0.082** | 0.500 |
| k562 | unseen (300) | mean+drop | 0.03 | 0.103 | 0.499 |

- Result: transfer is best on correlation (0.22–0.39) and discrimination (0.34–0.41) in every held-out line.
- mae: transfer is about equal to "no change" in rpe1/hepg2/jurkat, but clearly worse in h1 and k562.
  h1 is the most different line (stem cells, VCC measurement). 2026 uses yet another measurement (10x Flex).
- For genes never switched off elsewhere ("unseen"), nothing beats "no change" by much.
- The "seen" genes here are mostly essential genes with strong effects. The 2026 targets are not essential,
  so their changes are probably weaker; these numbers are likely optimistic for 2026.

### Exp 5: first 2026 leaderboard submission (transfer, shrunk by half)
- `src/predict_2026.py` → `data/submissions/transfer_a05.vcc`
- For each of the 300 targets: average its change over the public lines where it was switched off
  (272 in K562; 25 also in H1), shrink by 0.5, make the switched-off gene drop to ~0.31x.
  The 28 never-seen targets get the average change, shrunk the same way.
- Applied to each context's control cells; 400 new cells per target drawn around the prediction
  (totals and per-gene spread taken from that context's controls; no control cell copied).
- Shrink factor 0.5 from Exp 4: vs no shrink, error dropped in 4 of 5 held-out lines, correlation unchanged.

| held out | mae, no shrink | mae, shrink 0.5 | mae, no change |
|---|---|---|---|
| h1 | 0.719 | 0.496 | 0.406 |
| k562 | 0.158 | 0.126 | 0.133 |
| rpe1 | 0.273 | 0.263 | 0.281 |
| hepg2 | 0.334 | 0.312 | 0.338 |
| jurkat | 0.239 | 0.219 | 0.228 |

- Arc's `vcc prep` needs ~33 GB RAM for a full submission, so the script writes and packages the file
  in chunks itself. A 6-gene pilot passed `vcc prep --dry-run` and the .vcc archive check.
- Scoring: 6 metrics per context, 0 = Arc's mean-response baseline for that context. 2 submissions/day.
- File built and checked: 360,000 cells x 18,533 genes, 400 cells per target per context, whole counts,
  targets match the official list, ~5.4k nonzeros per cell (real controls ~6.0k), 3.0 GB .vcc.
- Submitted 2026-09-30 as `ask_sathvik` (entry wgm3kKK92ecNULbKTbtn). Upload 20 min, scoring a few min.

| overall | pds | mse | nmae | fid | reach | jac | rank |
|---|---|---|---|---|---|---|---|
| **0.097** | 0.506 | 0.000 | 0.070 | −0.088 | 0.124 | −0.031 | 544 |

- Scores are scaled: 0 = Arc's average-change baseline (built from the true answers), 1 = a replicate experiment.
- Result: above the baseline overall, almost all of it from telling genes apart (pds 0.51).
  Expression error (mse) is at its floor of 0; the direction of significant genes (fid, jac) is worse
  than the baseline.
- For reference, a public write-up of the same idea without shrinking scored 0.046 (pds 0.41, fid −0.18).

---

## Current

### Codex transfer-calibration experiments
- New calibration and sampler experiments build on Exp 5/6 without replacing the original predictor.
- Details, results, scoring-space correction, and commands: [CALIBRATION.md](CALIBRATION.md).
- HepG2 expression shrinkage 0.25 is rejected: official overall 0.009691 versus 0.170889 at 0.5.
- At shrinkage 0.5, correlated noise scores 0.182908 and depth-aware gamma-Poisson 0.184269 on HepG2, versus 0.170889 original. Gains need seed confirmation.
- Testing depth-aware noise on Claude's stronger alpha=1 baseline; a calibrated submission packager's small pilot passes Arc's dry-run.

### Exp 6: local copy of the official scorer
- Goal: try ideas without spending leaderboard submissions (2/day).
- `src/local_eval.py`: a held-out public line plays a 2026 context. 150 random targets, up to 200 cells each.
  Its control cells are split in half: one half is the model's input, the other goes with the answers.
- Scored with Arc's scorer (`cell-eval2`, preset `vcc2026`): official baseline = 0, replicate = 1.
- Lines: HepG2 done, H1 running (bigger: 18k genes, baseline alone takes >30 min).
- The prediction file must also hold control rows; the scorer ignores them and uses the real held-out controls.
- `src/noise.py`: per-gene noise of each line, for denoising (weight each copied change by z²/(z²+k)).

HepG2 (150 targets), official scale (0 = baseline, 1 = replicate):

| method | overall | pds | mse | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|---|
| transfer, shrink 0.5 (= leaderboard 0.097) | **0.171** | 0.83 | 0 | 0.00 | 0.23 | 0.05 | −0.08 |
| same, other sampling seed | 0.156 | 0.81 | 0 | −0.01 | 0.21 | 0.09 | −0.17 |
| + denoise k=4 | 0.109 | 0.86 | 0 | −0.27 | 0.05 | −0.09 | 0.10 |

- Local scorer runs; same pattern as the leaderboard (discrimination carries the score, mse 0).
- Sampling seed alone moves the overall by ~0.015.
- Denoising hurts on HepG2 (changes copied from strong-effect lines). Waiting on H1, the more 2026-like case.

Finding: for the 2026 targets, K562's profiles are mostly noise.
- Genes changed by more than 5 noise units, per knockdown (switched-off gene not counted), median [90th pct]:
  K562 all 1 [15], K562 2026 targets 1 [8]; RPE1 20 [268]; HepG2 4 [172]; Jurkat 5 [59];
  H1 309 [3049], H1 2026 targets 566 [3133] (H1 has ~1,000 cells per knockdown vs ~100–150 in K562,
  and its noise is per-cell only, without batch effects, so H1 counts are inflated).
- Common response (average over all knockdowns) is only partly shared between lines: corr 0.3–0.9
  (H1 vs others ~0.3–0.4). The average of the other lines predicts a line's common response no better
  than K562 alone (H1 0.38 vs 0.40, Jurkat 0.45 vs 0.61).
- Trying (branch `exp/calibration`): weight each line's copied change for a knockdown by s/(s+kt),
  s = number of clearly changed genes, so near-empty profiles shrink to ~0.

Shrinking vs not (local official scorer; HepG2 on the official scale, H1 relative to baseline only):

| method | HepG2 | H1 |
|---|---|---|
| transfer, shrink 0.5 (leaderboard 0.097) | 0.171 | 0.177 |
| + denoise per gene (k=4) | 0.109 | 0.130 |
| + weight per knockdown (kt=3) | 0.122 | — |
| no shrink (alpha 0.75) | 0.282 | |
| **no shrink (alpha 1)** | **0.374** | 0.178 |
| alpha 1.25 | 0.347 | |
| alpha 1.5 | 0.291 | |
| alpha 2 | 0.117 | 0.069 |

- Any shrinking hurts: mse is already at its floor of 0, so smaller changes gain nothing there and only
  lose significant, correctly-signed genes (fid, reach). Past alpha 1, nmae (error on the size of changes) collapses.
- Generated no-change cells vs real controls: 4 genes falsely called changed (real vs real: 0), so the cell
  generator is not the problem.

Leaderboard check of alpha 1 (entry 9Uk2yux4d7EyfOtcjN7Y, 2026-09-30):

| submission | overall | pds | mse | nmae | fid | reach | jac | rank |
|---|---|---|---|---|---|---|---|---|
| alpha 0.5 | 0.097 | 0.506 | 0 | 0.070 | −0.088 | 0.124 | −0.031 | 544 |
| alpha 1 | 0.092 | 0.509 | 0 | −0.048 | −0.026 | 0.137 | −0.018 | 556 |

- No gain on the leaderboard, as H1 predicted (0.178 vs 0.177), not HepG2 (0.374 vs 0.171).
  The components moved the way H1 showed: direction scores up, nmae down. H1 is the local test to trust.
- The leaderboard keeps only a team's newest submission, so it now shows 0.092.


Removing the part shared by all targets (gamma 0: subtract the panel-average change, keep the switched-off gene's drop), H1 relative to baseline:

| method | H1 | pds | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| alpha 0.5 | 0.177 | | | | | |
| alpha 0.5, gamma 0 | 0.173 | 0.64 | −0.04 | 0.35 | 0.04 | 0.05 |
| alpha 1 | 0.178 | 0.66 | −0.16 | 0.47 | 0.02 | 0.08 |
| alpha 1, gamma 0 | 0.176 | 0.66 | −0.16 | 0.45 | 0.02 | 0.08 |

- No gain; within seed noise (~0.015). The shared part neither helps nor hurts discrimination. Not submitted.

Sampler spread and seed noise on H1 (alpha 1, relative to baseline):

| method | H1 | pds | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| alpha 1 (seed 0) | 0.178 | 0.66 | −0.16 | 0.47 | 0.02 | 0.08 |
| alpha 1, seed 1 | 0.177 | 0.65 | −0.16 | 0.46 | 0.02 | 0.08 |
| alpha 1, no extra spread (phi 0, Poisson) | 0.179 | 0.65 | −0.15 | 0.48 | 0.01 | 0.09 |

- On H1 the seed moves the score by only ~0.001 (HepG2: ~0.015). Plain Poisson makes no real difference. Not submitted.

Finding: the public lines barely cover the 2026 targets. K562 covers 272 of 300, H1 25, RPE1/HepG2/Jurkat 0.
The 2026 contexts' control cells correlate with H1 at 0.50/0.63/0.71 (A/B/C) and with K562 at 0.33/0.22/0.31.
So the typical knockdown response we copy comes from the wrong cell type. Testing H1's typical response
plus K562's target-specific part (math in [MODEL.md](MODEL.md)).

H1 typical response + K562 specific part (MODEL.md), H1 relative to baseline:

| method | H1 | pds | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| alpha 1 (bar) | 0.178 | 0.66 | −0.16 | 0.47 | 0.02 | 0.08 |
| alpha 1 + H1 typical | 0.162 | 0.47 | −0.12 | 0.47 | 0.06 | 0.08 |
| alpha 0.5 + H1 typical | 0.134 | 0.33 | 0.01 | 0.38 | 0.03 | 0.05 |

- Worse. Adding the same typical change to every target makes their predictions alike, so telling them apart (pds) drops sharply.
  The significant-gene scores gain little (reach +0.04). pds wants only each target's own part.

Fast offline check (seconds, not an hour): 280 H1 targets that are also in K562. We predict H1's change from K562's,
and for each target count the share of other targets whose true profile is closer than its own
(cosine, log1p CP10k, own gene excluded; lower is better; this is close to pds):

| K562 profile used | share closer (lower better) |
|---|---|
| raw | 0.133 |
| centered (gamma 0) | 0.124 |
| averaged with its 5–50 most similar knockdowns | 0.157–0.161 |
| genes with abs(z) < 1 / 2 / 3 set to 0 | 0.132 / 0.142 / 0.170 |
| capped at ±0.5 / ±1 | 0.141 / 0.133 |

- Every cleanup that pulls profiles together or drops weak genes makes targets harder to tell apart.
  The raw (or centered) K562 profile is as good as it gets from K562 alone.

Neighboring genes drop too (genes starting within 1 kb of the target set to the typical own drop, within 5 kb to half of it;
`add_neighbours`, gene positions from UCSC refGene). H1, relative to baseline:

| method | H1 | pds | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| alpha 1 (bar) | 0.178 | 0.664 | −0.16 | 0.47 | 0.02 | 0.08 |
| alpha 1 + neighbors | 0.180 | 0.677 | −0.16 | 0.47 | 0.02 | 0.08 |
| alpha 1, gamma 0 + neighbors | 0.180 | 0.688 | −0.16 | 0.46 | 0.02 | 0.08 |

- pds up 0.01–0.02. Overall +0.002, about twice the seed noise (~0.001). Small but consistent; kept as a default add-on, not submitted alone.
- Fast-check intervals (136 targets, paired bootstrap): K562 → H1 linear map −0.008 [−0.026, +0.008]; neighbors 0.000 [−0.024, +0.026].
  The fast check is too noisy for changes this small.

Next source: X-Atlas/Orion (Xaira 2025), genome-wide CRISPRi in HCT116 and HEK293T, per-cell parquet on Hugging Face (~127 GB).
Each batch mixes all targets (one test batch: 18.5k cells, 10k targets, 168 of the 2026 targets).

X-Atlas/Orion added (`src/xatlas.py`, streamed 127 GB, kept per-target means): HCT116 18,293 targets (median 150 cells),
HEK293T 18,311 (median 200). Both cover all 300 2026 targets.
Fast check: 279 H1 targets that all three lines knocked down; centered profiles; share of other targets ranked
closer (lower is better); paired bootstrap 95% interval vs K562:

| source | share closer | vs K562 |
|---|---|---|
| K562 alone | 0.124 | |
| HCT116 alone | 0.222 | +0.098 [+0.063, +0.131] |
| HEK293T alone | 0.241 | +0.116 [+0.083, +0.148] |
| HCT116 + HEK293T | 0.191 | +0.067 [+0.038, +0.096] |
| K562 + HCT116 | 0.121 | −0.003 [−0.023, +0.016] |
| K562 + HEK293T | 0.130 | +0.006 [−0.010, +0.021] |
| K562 + HCT116 + HEK293T | 0.118 | −0.006 [−0.026, +0.014] |

- No gain. The new lines are much worse than K562 at predicting H1, and adding them to K562 is noise.
- HCT116 shows ~34 clearly changed genes per 2026 target (K562 ~1), but those changes are cell-type specific:
  K562 and HCT116 profiles of the same target barely agree (the right target ranks behind 26% of others).
- Knockdown effects are mostly cell-type specific, so more unrelated cell types don't help. Not run on the full H1 test.
- Caveat: the 2025 H1 targets may have been chosen with K562 in mind, which would favor K562 on this check.

Gradient boosting cross-line model (`src/gbm.py`): 2.5M rows (destination line, target, gene) over K562, HCT116, HEK293T,
RPE1, HepG2 and Jurkat; label = destination's centered change; 2026 and H1 targets excluded from training.
Fast check: cosine to H1's true profile 0.078 → 0.121; pds-like −0.010 [−0.036, +0.016]. Full H1 test:

| method | H1 | pds | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| alpha 1 + neighbors (bar) | 0.180 | 0.68 | −0.16 | 0.46 | 0.02 | 0.08 |
| GBM changes | 0.089 | 0.54 | −0.06 | 0.06 | −0.01 | 0.01 |
| GBM rescaled to K562 size | 0.162 | 0.63 | −0.15 | 0.38 | 0.05 | 0.06 |

- Worse. The learned changes are closer to the truth on average but blander: targets look more alike (pds)
  and fewer genes come out clearly changed in the right direction (fid). Same lesson as Codex's learned
  shrinkage / PCA-ridge. Dropped.

Public #82 recipe on our sources (`src/atlas_shift.py`; generator vendored from kaipengm2/Virtual-Cell-Challenge-2026, MIT):
two targets per knockdown, per-cell-normalized mean scaled 0.6 and pooled profile scaled 0.3 (clipped at ±3);
centered sources weighted K562 2, H1 2, HCT116 1, HEK293T 1 (H1 excluded locally); promoter-neighbor prior;
cells = 4-control-cell templates reweighted to hit both targets exactly. H1, relative to baseline:

| method | H1 | pds | mse | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|---|
| alpha 1 + neighbors (previous best) | 0.180 | 0.68 | 0 | −0.16 | 0.46 | 0.02 | 0.08 |
| vetted K562 controls + neighbors, alpha 1 | 0.178 | 0.66 | 0 | −0.16 | 0.47 | 0.02 | 0.08 |
| **#82 recipe (0.6 / 0.3)** | **0.205** | 0.71 | 0.009 | −0.02 | 0.44 | 0.02 | 0.07 |

- First clear gain: +0.025 (~25x seed noise). It comes from pds, nmae and, for the first time, mse above 0.
  Scaling the two targets separately avoids the shrink trade-off we kept hitting with a single alpha.
- Vetted controls + neighbor rule don't stack (0.178).

Tuning the two scales (H1, relative to baseline):

| per-cell scale (ac) | pooled scale (ab) | H1 | pds | mse | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|---|---|
| 0.6 | 0.3 | 0.205 | 0.709 | 0.009 | −0.015 | 0.436 | 0.022 | 0.069 |
| 0.6 | 0.5 | 0.208 | 0.726 | 0 | −0.016 | 0.440 | 0.025 | 0.071 |
| 1.0 | 0.5 | **0.217** | 0.727 | 0.001 | −0.024 | 0.470 | 0.044 | 0.087 |

- The per-cell scale drives the significant-gene scores (fid, reach, jac); 1.0 beats 0.6 by ~0.01.
  The pooled scale trades mse against pds. Testing ac 1.5.
- ac 1.5 fails: the generator can't hit both targets when they pull too far apart ("Moment fitting failed").
  ac 1.0 / ab 0.5 is the setting.
- Seed check of ac 1.0 / ab 0.5: seed 1 = 0.2185 vs seed 0 = 0.2174 (same parts within 0.007). Confirmed.
  Built as `data/submissions/as_c10_b05.vcc` (validated: 360k cells, 300 targets x 400 x A/B/C, 2.42e9 nonzeros).

Leaderboard, 2026-10-01: `as_c10_b05` (entry f16jrONNUjaUMWbq9HIl)

| submission | overall | pds | mse | nmae | fid | reach | jac | rank |
|---|---|---|---|---|---|---|---|---|
| alpha 0.5 (old) | 0.097 | 0.506 | 0 | 0.070 | −0.088 | 0.124 | −0.031 | 544 |
| alpha 1 (old) | 0.092 | 0.509 | 0 | −0.048 | −0.026 | 0.137 | −0.018 | 556 |
| **#82 recipe, ac 1.0 / ab 0.5** | **0.159** | 0.637 | 0.073 | 0.156 | −0.019 | 0.105 | 0.000 | **292** |

- +0.062 over the previous entry (0.092), +0.062 over our best (0.097). Above the public #82 entry (0.155),
  without their CD4 and per-cell K562 sources. H1 predicted the gain (0.180 → 0.217).
- Gains: pds +0.13, mse 0 → 0.07, nmae +0.09–0.20. fid and jac are still about 0 (no better than the baseline):
  the significant-gene sets are where the remaining room is.

Why fid/jac are ~0 (H1, current best ac 1.0 / ab 0.5, per-gene DE tables saved with --cache-pred):

| per target | median | mean |
|---|---|---|
| genes we make significant | 4,270 | 4,250 |
| genes truly significant | 268 | 1,239 |
| share of our significant genes with the right sign | 0.52 | 0.54 |

- We call ~16x too many genes, with near coin-flip signs, and our count doesn't track the true count (corr 0.01).
- Likely cause: each predicted cell averages 4 real control cells, so predicted cells have fewer zeros and less
  spread than real ones. The Wilcoxon rank test picks that up on thousands of genes even where the mean is right.
- Cutting our list to the true length wouldn't help (fid 0.517 → 0.533): the artifacts sit at the top of the ranking.
- Testing templates from single cells (pool 1) and pairs (pool 2).
- Pool size (H1): pool 4 = 0.217, pool 2 = 0.217, pool 1 = 0.206 (mse 0.045, reach 0.077, but fid 0.36, jac 0.05).
  The generator's reweighting bends every gene a little, and the direction of that artifact depends on pool size:
  share of our significant genes called UP is 81% / 66% / 17% for pool 4 / 2 / 1 (truth 48%), with sign precision
  0.535 / 0.544 / 0.550 and median count 4,270 / 2,041 / 842 (truth 268).
- Testing binomial thinning instead (Gerard 2020, seqgendiff): real control cells, each gene's counts scaled by the
  per-cell target's ratio to the control mean; unchanged genes stay exactly as the real cells. Runs `as_thin_c10`, `as_thin_c06`.
- Binomial thinning: worse (0.195 at ac 1.0, 0.157 at ac 0.6; pds drops to 0.65 / 0.59). It still calls a median of 1,764 genes
  with 0.539 sign precision, so the template artifact is not the main problem: our transferred changes themselves
  have the wrong sign in H1 about half the time.
- Offline (H1, 150 targets): sign accuracy of transferred changes vs truth, on truly significant genes:
  all 0.61; |change| > 0.25 and all 3 sources agree 0.76 (~236 genes per target, close to the true median 268);
  |change| > 0.5 and 3 agree 0.82 (22 per target).
- Testing: per-cell target keeps only changes > 0.25 where all 3 sources agree; pooled target unchanged
  (`as_ag3_p4`, `as_ag3_p2`).
- Source-agreement mask on the per-cell target: worse (0.161 pool 4, 0.149 pool 2; fid 0.26 / 0.17 vs 0.47).
  fid also rewards how many correct calls we make: fewer calls gave fewer correct ones while the true count
  still sets the denominator. Keeping every change (current best) stays best.
- Codex: CD4 T-cell source (weight 0.5, `src/atlas_cd4.py`) at ac 1.0 / ab 0.5: H1 0.2235 vs 0.2174 / 0.2185 (two seeds)
  = +0.005; pds +0.017, mse +0.015, fid/reach/jac flat.
- Gene-similarity idea (predict change only for genes similar to the target): offline AUC for picking H1's truly
  significant genes per target (132 targets): expression level 0.646; transferred change x sqrt(expression) 0.607;
  co-expression with target in input controls 0.549; size of transferred change 0.529; K562 knockdown-profile
  similarity 0.507; general responsiveness in sources 0.481. Which genes become significant is mostly statistical power
  (expression), which our generator already gets from real control cells. Similarity adds little. Not pursued.

Overnight 2026-10-01 (chain v2 + Codex's shared selector `src/overnight_select.py`):
- K562 from the per-cell file (batch-matched controls, `src/k562_cells.py`): H1 0.219 vs 0.217 / 0.219 for the averaged
  file, a tie. The pooled-share stand-in was good enough.
- K562 per-cell + CD4 (w 0.5): 0.2253. Codex's CD4 at weight 1 (averaged K562): 0.2269 / 0.2276 (two seeds), the best.
- Selected and uploaded `as_cd4_w1_c10_b05` (entry KHbUyq5agxuOh2mSHREg):

| submission | overall | pds | mse | nmae | fid | reach | jac | rank |
|---|---|---|---|---|---|---|---|---|
| as_c10_b05 | 0.159 | 0.637 | 0.073 | 0.156 | −0.019 | 0.105 | 0.000 | 292 |
| **as_cd4_w1_c10_b05 (+ CD4 w 1)** | **0.1745** | 0.667 | 0.103 | 0.174 | −0.012 | 0.115 | 0.000 | **232** |

- +0.016 on the leaderboard from adding CD4 (H1 had +0.009). fid/jac still ~0.

### Codex functional perturbation graph — October 2

Literature-backed source-only cosine graph transport implemented in src/perturbation_graph.py; math and sources in AMBITIOUS_BETS.md. Reference set 747 shared cached training knockdowns, all H1/2026 targets excluded. Held-out summary screen on H1 (279 targets), HepG2/Jurkat (98 each), positive cosine neighbors k5/20/50, blends .25/.5/1, shrink and norm-matched controls. Rejected: k20/blend.5 H1 MSE .006652 versus half-copy .006546, retrieval .867088 versus copy .882812; same error/discrimination failure in other two lines. Norm-matched H1 MSE .010164 vs copy .010249 is tiny and retrieval still declines. These are common-gene centered-logratio proxies, not official metrics, not a comparison against full CD4 baseline. Two transport tests pass. No full scorer or upload. Saved data/calibration/perturbation_graph_screen.csv.

## 2026-10-02: per-context source weights (Jurkat vs HepG2)
- A/B/C markers: A T-cell-like, B mesenchymal, C squamous epithelial; none stem-cell.
- New local eval `jurkat`. CD4 weight 0/1/4 → Jurkat .1946/.1976/.1938, HepG2 .2365/.2314/.2165.
- Cell-type matching is real but small (±.003–.005). Target's own expression doesn't predict response size (dead end).

### Codex iPSC expression source — October 2 (full results pending)

Feng/HipSci figshare26819743 GenomeWideScreen_LFC_byGene.tsv.gz downloaded (794610300bytes), official MD5 verified. 43187656rows/6673targets; 6151 measured challenge genes. Source LFC is log10, converted by log2(10); adjusted mean log-expression effects, not exact pooled-count ratios. Streaming source src/ipsc_source.py writes ipsc_de.npz and separate ipsc_de_eval.npz (extra HepG2/Jurkat targets). Main challengecoverage182/300 all,77 ownQC; H1coverage250/300 all,134QC. Full-panel centering excludes own coordinates, missingness retained, ambiguous symbols/duplicate mapped pairs rejected. Five tests pass, including state restoration; corrected screen weight0 matches CD4w1 exactly on H1/HepG2. Initial limited-target HepG2 proxy preserved and superseded by expanded atlas_shift_x proxy. Conservative full candidate ownQC weight.5/ac1/ab.5/pool4: H1 and HepG2 as_ipsc_qc_w05_c10_b05 running, max2total; matched baselines .2268683314(H1) and .2313675357(HepG2 cw_base). Math/provenance/commands: IPSC_SOURCE_MATH.md. Claude owns uploads; none by Codex.


## Completed iPSC decision — October 2

| Local comparison | Baseline | iPSC w.5 | Difference |
|---|---:|---:|---:|
| H1 seed0 pilot | .226868 | .236112 | +.009244 |
| H1 corrected seed1 | .227553 | .237081 | +.009529 |
| HepG2 corrected seed0 vs CD4w1 | .231368 | .234210 | +.002843 |
| Jurkat corrected seed0 vs CD4w1 | .197603 | .192779 | −.004823 |

HepG2's tested no-CD4 model scores .236471, above this iPSC blend. All six Jurkat metrics decline. H1 gains persist when excluding panel target genes, but do not establish transfer to non-stem contexts. Reject the uniform iPSC addition for submission. Seed0 H1 is the pilot implementation; seed1 and the corrected other-line runs use exact numerical isolation. Compare scores within each line because normalization differs. All Codex scorer jobs finished, both reservations released, and Claude remains sole upload owner. No iPSC submission was built or uploaded. A no-CD4+iPSC ablation remains untested. Six source/adapter tests passed; cache and commands above remain reusable. This completed decision supersedes earlier pending-run snapshots.

## 2026-10-02: leaderboard probe ctxC_hct → 0.1706 (vs 0.1745)
- Same as as_cd4_w1_c10_b05 except context C: HCT116 ×2, no CD4. Only C changed, so C's own score fell ~.012.
- pds .655 (was .667), mse .099, nmae .166, fid −.014, reach .119, jac .000.
- Lesson: CD4 helps the epithelial context too. CD4's gain is not T-cell-specific, and HepG2 is a poor proxy for C.
  Pooling diverse sources beats cell-type matching here. KOLF (stem) still only helps H1 (paired r), not added.


2026-10-02 update: full VIPerturb blend scored 0.193075 on Jurkat versus 0.197603 baseline and is rejected for submission. Bulk-only ratio 0.25 ablation is now running (as_vip_bulk_r025_c10_b05); all 150 desired per-cell profiles were verified bitwise equal to baseline before launch. This isolates pooled expression changes. User owns submission approval.


2026-10-02 final bulk-only VIP Jurkat result: 0.1924778901 versus cw_base 0.1976025272 (delta -0.0051246371). Normalized MSE 0.082936 vs 0.091138; REACH 0.081389 vs 0.102610. Proxy improvement did not translate. Reject uniform full and bulk-only VIP blends for submission. Scorer complete, slot released.

- 2026-10-02 tv125 FINAL Jurkat .1869454 vs .1976025 baseline (reject). FID .107894 vs .212005; REACH .135591 vs .102610; JAC .030032 vs .020601; MSE .091141 vs .091138; PDS .702225 vs .702317; NMAE .054790 vs .056945. Strong distribution tradeoff. Previous scorer exited0; preflight shows no live scorer. Reserving ONE slot for opposite-direction tv075 (power.75), same profiles/depths.

## 2026-10-02: agreement allocation → leaderboard **0.1845** (rank 208), from 0.1745
- aa75 = as_cd4_w1_c10_b05 + per-target bulk amplitude × A^.75 m^-.75 (A = mean pairwise source cosine), energy kept,
  cap 2× (`src/agree_alloc.py`).
- LB: pds .669 (=), **mse .159 (from .103)**, nmae .174, fid −.013, reach .118, jac .001.
- Local: H1 +.002, Jurkat +.002, HepG2 −.012. Local tests under-weight mse; the leaderboard rewards it.

- 2026-10-02 CORRECTION + RESULT: tv075 finished naturally before attempted SIGSTOP; kill returned no such process, so it was NEVER paused. Session exited0. Jurkat tv075 .2149141711 vs cw_base .1976025272 (+.01731164). FID .302616 vs .212005; REACH .113110 vs .102610; JAC .027337 vs .020601; MSE .091099 vs .091138; PDS .702591 vs .702317; NMAE .052730 vs .056945. Promising distribution lever, needs HepG2 and seed replication / aa75 combination before proposing upload. No Codex scorer live; your two jobs may still occupy capacity.

## 2026-10-02 afternoon: amplitude, variance, source weighting (all on top of aa75)
- Bulk amplitude ab: H1 .35/.5/.75/1.25 → .2267/**.2288**/.2230/.2211. The Jurkat gains at large ab (.2750 at 1.25)
  are a scorer artifact: the sampling-correction deduction grows with spread, and plain mse got worse on 82% of
  targets. Keep ab .5.
- **Per-cell amplitude ac 1.25** (fail-soft generator): H1 **.2365**, Jurkat **.2192** (aa75: .2288/.1996). Plain mse
  better on 85% of targets. ac 1.5 ties on H1.
- Dropped:
  - per-cell agreement scaling .2253
  - pool 2: .2335 / .1937
  - per-target source agreement weighting γ1/γ2: .2250/.2121
  - Codex tpow .75 + aa75 on H1: .2269 (cells underdispersed)
  - α .5 / β 0: tie
- Evening, on top of aa75_c125 (H1 .2365 / Jurkat .2192), all ties or losses:
  - Flex K562 (VIP union bulk-only r.25): .2378 / .2185
  - per-cell norm restoration nr .5: .2378 (H1); nr 1 @ ac 1: .2354
  - per-gene sign-agreement reweighting gconf 1/2: H1 .2339/.2311; Jurkat gconf 1 .2175 (reach .153 → .124)
- Leaderboard raw metrics: the gap to #100 (.215) is reach (.184 vs .261 raw) and pds (.802 vs .824).
  fid/jac are hard for everyone; our mse matches #100.
- Scorer reach = among the REAL significant genes, ranked by our |predicted LFC| (our significant ones first), the
  deepest ≥90%-sign-pure prefix / N_conf. So it needs correct signs on our largest predicted changes.
- Night, on top of aa75_c125:
  - per-cell shrink of weakly expressed genes (ew 500): **H1 .2003 / Jurkat .1700**. Offline re-ranking looked good,
    but shrinking kills significance calls.
  - boost of well-expressed genes ×(1+.5√(cpm/(cpm+500))): H1 .2364 / Jurkat .2311; boost 1: H1 .2340.
  - **ac 1.5: H1 .2365 (tie) / Jurkat .2319**, so it's the 2nd-slot candidate. Built aa75_c15.vcc (0 fallbacks).
- Offline reach simulator `src/reach_sim.py` (matches the scorer within ~.003 on H1).
- Late night, on top of aa75_c125 (H1 .2365 / Jurkat .2192):
  - panel centering of the pooled change pc .5 / 1: .2371 / .2350 (H1), so a tie
  - **KOLF w1: H1 .2591 / Jurkat .2180**. Stem→stem gain, tie on Jurkat; candidate for a later day (55/300 panel).
  - X-Atlas HCT116/HEK293T weight 2: .2275 / .2019, worse.

## 2026-10-03 leaderboard
| submission | overall | rank | pds | mse | nmae | fid | reach | jac |
|---|---|---|---|---|---|---|---|---|
| aa75_c125 (ac 1.25) | **.1877** | 201 | .671 | .162 | .168 | −.010 | .133 | .002 |
| aa75_c15 (ac 1.5) | .1860 | 204 | .676 | .166 | .141 | −.009 | .139 | .002 |

ac ≈ 1.25 is optimal on the leaderboard: bigger ac helps pds/reach/mse but overshoots fold changes (nmae). The board
responds much more weakly than local tests (+.003 vs H1 +.008 / Jurkat +.02).
- 2026-10-03 morning: KOLF + VIP bulk r.25 on aa75_c125: H1 .2586 / Jurkat .2143 vs KOLF alone .2591/.2180. VIP adds
  nothing. kolf_c125.vcc (validated) is scheduled for 2026-10-04 05:32 IST.
- 2026-10-03: DE-consensus boost/suppress on the per-cell change (top-10% agreeing genes ×1.6, rest ×.5, per #296's
  description): H1 .2252 / Jurkat .1961, worse (fid and reach down). Per-target H1 pds median .98 (retrieval easy on H1),
  while LB raw pds is .80, so H1 can't guide pds work.
- 2026-10-03: agreement alpha 1.5 (capped) + ac 1.25: H1 .2352 / Jurkat .1925 (pds .552), worse. Keep alpha .75.
- 2026-10-03 **experimental stem recipe** (October 4 correction: Arc says the 2026 lines differ from 2025 H1; the blog does not establish a stem context in the final panel):

  | H1 | overall | pds | mse |
  |---|---|---|---|
  | aa75_c125 | .2365 | .747 | .078 |
  | +KOLF | .2591 | .787 | .127 |
  | +iPSC (Codex's HipSci 34-line, own-suppression QC) | .2494 | .799 | .093 |
  | **+KOLF+iPSC** | **.2673** | .814 | .139 |

  `src/final_router.py` scores pluripotency markers in each context's controls: H1 6.54 → stem; A/B/C .10–.39 and
  K562/Jurkat/HepG2 0 → default. Testing KOLF 2 + iPSC 2.
- 2026-10-03: stem weights KOLF/iPSC 1/1, 2/2, 4/4 → H1 .2673 / **.2712** / .2672; router stem recipe set to 2/2.
- 2026-10-03: iPSC w1 on Jurkat (non-stem): .2131 vs .2192, worse (KOLF tie .2180). Stem sources only for stem contexts, which supports the router split.
- 2026-10-03: default recipe without the H1 source, on Jurkat: .2154 vs .2192. Keep H1 weight 2.

## 2026-10-04 Codex neural transfer and external-feature experiments

Actual aa75_c125 pooled predictions form the baseline for a rank32 response residual network. Outer held line is excluded from training labels and source inputs; RPE1 selects checkpoints. Model development has used aggregate held-line results, so these are development tests, not untouched final estimates.

| Full generated-cell comparison | Baseline | Candidate | Interpretation |
|---|---:|---:|---|
| Jurkat, uncentered residual | .219159 | .245119 | PDS improves; uncapped unbiased expression MSE worsens |
| Jurkat, centered residual | .219159 | .249920 | Raw PDS .83396→.87096; reach falls; uncapped MSE improves only36% targets |
| HepG2, centered residual | .263551 | .261908 | Matched replicate anchors; slight net loss |

HepG2 averages above use from_replicate; Jurkat uses from_baseline. Never compare their absolute magnitudes as a shared scale. HepG2 also loses on the mean of its six baseline-normalized metrics (.151770→.150169). Raw audits and normalization audit are in data/calibration/atlas_nn_*raw_audit.csv and neural_normalization_audit.json. No neural submission made.

Cheap controlled follow-ups, all completed:

- Within-context contrastive loss: mixed retrieval, slightly worse MSE than original centered model on both lines. Not promoted.
- Rank128: worse MSE and retrieval on both lines. Oracle centered-basis error coverage rises modestly, but learned transfer does not exploit it. Not promoted.
- ESM2 protein features: verified595,319,311-byte Altos/Arc source,300/300 challenge coverage. Training-only16-component PCA, explicit missing flag. Real vs identity-shuffled features essentially tie on Jurkat; HepG2 real improves retrieval but loses MSE. No consistent protein-specific benefit. See PROTEIN_FEATURE_BET.md and data/calibration/protein_feature_comparison.json.
- Jiang pathway source: archive checksum verified;1,626 line/target/pathway profiles,218 distinct targets, only9 challenge overlaps. Stimulated context and missing LFCs prevent treating this as a dense steady-state source. Diagnostic agreement with centered K562 source is weak; no blend promoted. See data/jiang/inventory_summary.json and transfer_screen.csv.

H1 frozen rank32 centered/no-protein check completed: full-profile MSE .00433089→.00435327, retrieval .858837→.860626, cosine .080757→.085665. Weak mixed evidence; no full scorer promoted.

## 2026-10-04 official KOLF probe

`submit_kolf_c125.log` confirms published entry EbA8m4YGpsYAR4SLWSYB: **.18594035**, rank213 at scoring, below aa75_c125 .18765290. The strong H1 local KOLF gain did not transfer to validation A/B/C. Best and latest submission differ. Codex made no upload.


### 2026-10-04 Confidence weighting: full generated-cell rejection

Fixed cross-context sign classifier trained on H1/HepG2, tested on Jurkat. Direct positive confidence scaling scored0.208691 vs aa75_c1250.219159. Restoring per-target nonself log-effect norm before count normalization recovered0.218415, still below baseline. Reach improved slightly (.152958→.157970), but NMAE/Jaccard deteriorated;0generator fallbacks. No submission and no further blind strength sweep. Cached hypothetical rank gains were not deployable gains. See data/calibration/aa75_signconf_s05_norm_jurkat.log and src/sign_confidence_diagnose.py.

## 2026-10-04 leaderboard
- kolf_c125 (aa75_c125 + KOLF w1): **.1859** (rank 213) vs .1877. pds .674, mse .154, nmae .167, reach .131. That's a tie,
  as Jurkat predicted: stem sources don't help non-stem A/B/C. Best stays aa75_c125 .1877.

## 2026-10-05 context PCA boost (new lever: destination control covariance)
- Real target-specific effects are concentrated in the destination controls' top PCs (100 PCs: 31% H1 / 18% Jurkat of
  the energy vs 7% / 13% in ours).
- Plain boost E + 4P over PCs 1–100 raises corr(pred, real) by 28%, but collapses retrieval:
  H1 full scorer .2197 (pds .747→.609).
- **Target-specific boost** (P from E minus the panel mean, PCs 4–100, g 2, norm kept; `pcs=100 pcg=2 pclo=3 pcspec=1
  pcn=1`):

  | | aa75_c125 | PCA boost | pds | mse |
  |---|---|---|---|---|
  | Jurkat | .2192 | **.2398** | .660→.732 | .124→.174 |
  | H1 | .2365 | .2333 | .747→.705 | .078→.102 |

  Plain mse better on 86% (Jurkat) / 63% (H1) of targets, so it's real. HepG2 check running.
- **Leaderboard: pca_c125 = .1765** (rank 269) vs .1877. pds .671→**.584**, mse .162→.183. HepG2 confirmed the Jurkat
  gain (.2636→.2710), but the **leaderboard followed H1** (pds down).
  - Lesson: for changes that move pds, trust H1 (same 10x Flex platform as 2026), not Jurkat/HepG2 (10x 3').
  - Earlier lessons still hold: for source choice, don't trust H1 (stem biology); for mse, check plain mse.
  - g 4 version not submitted.


### 2026-10-05 Distribution post-training (PerturbCellRL-inspired)

Implemented differentiable energy-distance endpoint fitting, then an ablation adding a pretrained endpoint anchor and a training-cell-calibrated expression-support tail penalty.200updates each; unchanged initial checkpoint included in validation selection. This adapts paper ideas to our deterministic transport and is not stochastic NFT. Real HepG2 perturbed cells withheld from fitting, identical test samples verified against original run, training-only representation replay confirmed.

Held150target program-space energy (lowerbetter): source-transfer0.821651; originalflow0.861163; energy-only1.186985; guarded1.102670. Mean and variance errors worsened too. Both seen-target93 and unseen-target57subgroups failed. Within-context held-target validation improved, showing that this training success did not transfer to a new cell context. Rejected for full count/scorer/upload promotion.7tests passed. See PERTURBCELLRL_EXPERIMENT.md and data/context_flow/hepg2_s0/posttraining_comparison.json.


### 2026-10-05 HepG2 context-PCA requested replication

Exact requested settings: alpha=.75 beta=.75 ac=1.25 pcs=100 pcl=.8 pcn=1 x=1. After iCloud stalls, run completed with0generatorfallbacks. Same normalization references verified. Replicate score0.263551→0.245268, baseline score0.151770→0.139086. Raw uncapped MSE improved26.9%, but PDS and reach degraded enough to dominate. Rejected for upload. Comparison data/calibration/pca_b_hepg2_comparison.json.

## 2026-10-05: remove shared control-PC directions (rmpc5) and more cell variance (tvar125)
- rmpc5 (`pcs=100 pcg=-1 pchi=5 pcn=1`: subtract each prediction's projection on the context's top-5 control PCs, renormalize): H1 .2408 vs .2365 (pds .747 -> .778; plain mse flat, .00290 -> .00292). pds gain concentrated on ~1/3 of targets.
- tvar125 (template power 1.25, from PerturbCellRL's under-dispersion point): H1 .2374, but the gain is capped-mse (.078 -> .090, gameable), fid/jac down. Dropped.

### Shared gene response rule (Codex, 2026-10-05)

- 2026-10-05 (Codex): New shared_response_rule.py tests a 10-feature, gene/target-ID-free modulation of copied effects by destination gene expression, target expression, and source agreement. Whole RPE1 context selects ridge; each outer context excluded. H1 desired-profile MSE .00500274 -> .00507881 (+1.52% worse), cosine .12309 -> .13070, retrieval slightly worse. Jurkat and HepG2 select zero correction. Rejected for promotion. No scorer/upload used. This and the centered-flow failure argue that simple expression-conditioned transfer is insufficient with current contexts.

Commands: `.venv/bin/python src/shared_response_rule.py --held {h1,jurkat,hepg2}` (run each context separately). Artifacts: `data/atlas_residual/<held>/shared_response_rule/`. Features all vanish when copied effects are zero; no generic additive response. Controls standardized with training-context statistics. Six ridge penalties selected on RPE1 only, with baseline as explicit candidate. Results are pooled-profile proxies, not emitted cells or official scores.

### Control covariance intervention pilot (2026-10-05)

- 2026-10-05 (Codex): Control-only covariance intervention tested across H1/Jurkat/HepG2. Normal-cell covariance supplies knockdown direction, matched to copied effect norm; blend/PC removal selected on another context. MSE improves 1.55% H1 and .26% HepG2, Jurkat slightly loses; retrieval loses all three. Norm-matched shrink audit: {"h1": {"baseline_mse": 0.005002736113965511, "network_mse": 0.00492521608248353, "norm_matched_shrink_mse": 0.004920904524624348, "direction_mse_delta": 4.3113027459185105e-06, "direction_delta_bootstrap95": [-1.0761611918042035e-06, 9.682663176135968e-06]}, "jurkat": {"baseline_mse": 0.02243981510400772, "network_mse": 0.022449087351560593, "norm_matched_shrink_mse": 0.022465776652097702, "direction_mse_delta": -1.668965160206426e-05, "direction_delta_bootstrap95": [-7.146589650801616e-05, 3.757478889383492e-05]}, "hepg2": {"baseline_mse": 0.03076029382646084, "network_mse": 0.03067927062511444, "norm_matched_shrink_mse": 0.03076963499188423, "direction_mse_delta": -9.036353003466502e-05, "direction_delta_bootstrap95": [-0.00012236298534844535, -6.275230098253817e-05]}}. No full scorer/upload. src/control_response_network.py; data/calibration/control_network_direction_audit.json. Inspired by control-only virtual knockout (https://pubmed.ncbi.nlm.nih.gov/35510185/), not scTenifoldKnk reproduction.

### Target-output protein geometry (2026-10-05)

- 2026-10-05 (Codex): Target-output protein geometry pilot completed (3 held contexts x real/shuffled annotation geometry). External ESM PCA16 supplies output-gene basis and target coordinates; ridge learns dynamics using source effect projections plus target features, generic-context correction centered, RPE1 selects penalty. Real versus shuffled MSE deltas: {"h1": {"real": {"mse_delta": 1.0800891363151738e-07, "retrieval_delta": 0.0004474272930649059, "selected_ridge": 1}, "shuffled": {"mse_delta": 1.6482746056506459e-06, "retrieval_delta": 0.0006263982102908683, "selected_ridge": 0.01}}, "jurkat": {"real": {"mse_delta": -1.427927194735945e-06, "retrieval_delta": 0.0013422818791946067, "selected_ridge": 0.1}, "shuffled": {"mse_delta": 2.8869886198601424e-06, "retrieval_delta": 0.0003131991051452676, "selected_ridge": 0.001}}, "hepg2": {"real": {"mse_delta": -7.088992027073915e-07, "retrieval_delta": 0.00031319910514537863, "selected_ridge": 1}, "shuffled": {"mse_delta": 4.3010457348446884e-09, "retrieval_delta": 0.0, "selected_ridge": 100}}}. Tiny differences only; no material candidate, no scorer/upload. src/protein_response_geometry.py; data/calibration/protein_geometry_comparison.json. Motivated by structure/dynamics separation in https://arxiv.org/html/2608.06824v1, not a GeneGeoFlow replication (their holdout is intervention, ours context).
- Jurkat (x=1; a first run without x=1 used the wrong caches and scored .001, discard): rmpc5 .210 vs .219 (pds .660 -> .705, but plain mse worse on 68% of targets, reach .153 -> .117). Offline proxy on Jurkat also improved retrieval (.848 -> .884); top-5 PCs hold only ~6% of prediction energy.
- rmpc5h (pcg=-.5): H1 .240 (pds .770), Jurkat .212. Same trade-off, smaller.
- Built and validated data/submissions/rmpc5.vcc (full strength). Open question it would answer on the board: does 2026 follow H1 (pds up, nothing else moves) or Jurkat (reach/mse down)?
- Reach diagnostic (aa75_c125, scratch reach_diag): sign accuracy in our top-10 ranked real-DE genes is .69 (H1) / .83
  (Jurkat), overall .62 / .73, vs .57 / .58 for "all down". k*=0 on 34/150 H1 targets. Reach is limited by sign accuracy
  at the top of our ranking, not by depth.
- 2026-10-05: uploaded rmpc5 (user OK). Expect ~.193 if the board follows H1, ~.181 if it follows Jurkat.

### Learned source-sign confidence (2026-10-05, Codex)

Gradient-boosted sign-correctness classifier uses copied effect/sign/magnitude, destination gene expression, target expression, source agreement/amplitude, and target-output protein cosine. Gene/target IDs excluded. Training samples up to128 gene effects per target with |observed pooled response| >= .1, contexts balanced; RPE1 chooses confidence exponent from 0/.5/1/2/4. All three choose4. Norm-preserved effect weighting: e_new = e*p(correct)^4, then restore original target effect norm.

- H1: top10 sign proxy .72378 -> .71911; MSE .00500274 -> .00491741; cosine .12309 -> .16271; retrieval .86416 -> .85839. Mixed.
- Jurkat: sign proxy .778 -> .79333; MSE .02243982 -> .02245210; retrieval .85951 -> .87396. Mixed.
- HepG2: sign proxy .77933 -> .80533; MSE .03076029 -> .03062613; cosine .15669 -> .16897; retrieval .86846 -> .87172. Promising.

These are desired-profile and large-response ranking proxies, not official DE reach. Full emitted-cell tests started on HepG2, names learned_sign_bulk (bulk only) and learned_sign_dual (bulk+per-cell). Both use fixed RPE1-selected exponent; no held-label tuning. Logs data/calibration/learned_sign_hepg2_{bulk,dual}_score.log. Full scorer results pending. No upload.

Early generated-cell mean-error audit for learned confidence: pooled-only .03268616 -> .03266475; dual .03268616 -> .03270225. Both paired-target bootstrap intervals include zero. This is plain pseudobulk MSE, no capped sampling correction. Full six-metric scores pending. Artifact `data/calibration/learned_sign_hepg2_generated_mse.json`.

### Confidence profile integration audit

- 2026-10-05 (Codex): Learned-sign integration audit: HepG2 global profile normalization reduces the common-axis MSE gain to ~2.0e-5 across native genes, and emitted pooled-only retains it (~2.1e-5); dual weighting reverses it. Generated desired-change cosine bulk .707/dual .737. Tried preserving modeled-subset mass as separate desired-profile diagnostic: retrieval improves all contexts vs global normalization, but MSE worsens on Jurkat/HepG2; H1 all metrics improve over baseline. Not a universal fix, no third scorer launched. data/calibration/learned_sign_hepg2_realization.json and learned_sign_mass_preservation.json. Existing bulk/dual scorer children still in progress; no upload.

### Learned-confidence full HepG2 results

- 2026-10-05 (Codex): Learned-sign HepG2 FULL results COMPLETE: dual improves replicate .26355075 -> .26962382 (+.0060731), baseline .15176974 -> .16884866 (+.0170789), identical verified anchor/bundle digests. Replicate metric deltas FID+.15538/REACH+.03272/PDS+.04111, NMAE-.08853/JAC-.10425; MSEscore0. One generator fallback at .8. Bulk-only mixed (replicate-.00321,baseline+.00250). Plain dual MSE slightly worse, CI includes0; gain is direction/discrimination scores, not capped-MSE artifact. Both completed prediction h5ads deleted by authorized cleanup. Taking next TWO free slots for unchanged dual H1 and Jurkat replication; no upload.

### Separate response-presence head (Codex)

- 2026-10-05 (Codex): Response-presence head tested three contexts, RPE1 selects power0 everywhere (leave sign-only rule unchanged). Presence is |true pooled delta|>=.1, not DE significance. Held AUC H1 .6943/Jurkat .6015/HepG2 .5987, poor context calibration: predicted presence H1 .229 vs real .083; HepG2 .359 vs .478. Rejected as weighting feature. src/response_presence_head.py; data/calibration/response_presence_comparison.json. Jurkat full dual sign run completed .219159 -> .224867, legacy anchors now being verified via fresh same-reference scoring. Jurkat slot RELEASED; H1 slot remains occupied. No upload.
- **rmpc5 leaderboard: .1943 (rank 223), new best** vs aa75_c125 .1877. pds .671 -> .713; mse .161, nmae .167, fid -.011,
  reach .134, jac .002 all flat. The board followed H1 again (pds up, nothing else moves), not Jurkat. H1 is the
  trusted predictor for pds-type changes on 2026. New default recipe: aa75_c125 + pcs=100 pcg=-1 pchi=5 pcn=1.
- **u3b15** = rmpc5 + boost ×1.5 the per-cell change of genes whose ≥3 sources all agree in sign (`cq=1 cboost=1.5
  csupp=1 cmin=3`; new `cmin` option). Motivation: sign accuracy among our top-50 ranked real-DE genes is .82 (H1) / .90
  (Jurkat) when ≥3 sources are unanimous vs .55-.70 otherwise; offline unanimous-first lifts Jurkat raw reach .308 -> .372.
  Unlike cons10 (Oct 3), nothing is suppressed.
  - Jurkat **.241** vs rmpc5 .210 / aa75_c125 .219 (fid .272 -> .372, reach .117 -> .187, nmae .076 -> .099).
  - H1 .242 vs rmpc5 .241 (tie; reach .053 -> .058, nmae -.032 -> -.051).
  - Plain mse better than rmpc5 on 99% (H1) / 97% (Jurkat) of targets.
  - HepG2: u3b15 vs aa75_c125: fid .216 -> .283, nmae .038 -> .058, pds .624 -> .628, reach flat (mean of 6 ~.152 -> .166). Built + validated data/submissions/u3b15.vcc (0 fallbacks); candidate for 2026-10-06 05:30 IST slot.

### 2026-10-05 learned confidence + neural combination

Full emitted-cell scores: H1 .242507, Jurkat .245837, zero generator fallbacks. Both exceed neural-half and rmpc5 locally; fresh scoring used identical SHA-checked baseline reference files. Comparisons: data/calibration/learned_sign_combo_shared_reference.json. Plain generated-cell MSE: data/calibration/learned_sign_combo_plain_mse.json. HepG2 replication completed: .270404 replicate / .159780 baseline versus .263551 / .151770 transfer baseline. Matching anchor/bundle digests, zero fallbacks. Plain MSE improves versus transfer baseline (paired interval excludes0); near-tie versus neuralhalf. No upload. Three-class signed response cheap pilot completed on all3: MSE improves but retrieval worsens throughout; not promoted. data/calibration/signed_response_comparison.json.
- u3b2 (cboost 2): Jurkat .256 (fid .450, reach .209, jac -.003) vs u3b15 .241; H1 .239 vs .242 (nmae -.051 -> -.077). Not a clear win on both; u3b15 stays first for 2026-10-06 05:31 IST. u3b2 built as 2nd-slot candidate if u3b15 lifts board fid/reach.

### 2026-10-05 source-aware signed-response pilot

Per-output-gene source effects/signs/missingness added for K562/H1/HCT116/HEK293T/CD4. Held context and each training example own source excluded. H1 MSE .005003->.004921, top10 large-response sign proxy .72378->.73178, retrieval .86416->.84957. RPE1-selected blend1. Not promoted due discrimination loss. Testing centered correction with same frozen classifier; no full scorer or upload. src/signed_response_head.py --source-features; artifacts signed_response_sources and signed_response_sources_centered.

### Source-aware signed head rejection and gated-flip hypothesis

Jurkat centered head worsens retrieval .85951->.81329 and MSE .022440->.022475. On substantial real responses, flipped signs are correct only .459 H1 / .486 Jurkat. Both raw and centered heads rejected for full scoring. Next cheap pilot keeps copied sign unless classifier conditional direction confidence exceeds RPE1-selected threshold .75/.9/1. Existing frozen classifier reused; held outcomes only reported. H1 gated pilot live. No upload/scorer added.
- Called-set diagnostic (u3b15, scratch callset.py), prompted by Codex's "Vcc2026 Ambitious Bet" note (layer 2):
  H1 calls median **4,644** genes per target vs 268 real (84% of targets over-call >1.5x), sign acc of calls .54;
  Jurkat 70 vs 66. Trimming our own calls by p-value to an oracle 1.2*n_real only lifts raw fid/jac .529/.094 ->
  .546/.112 (H1), .439/.066 -> .500/.158 (Jurkat). Test: `topk=K` keeps the per-cell change on each target's K largest
  |ec| (after the consensus boost), zero elsewhere, bulk change untouched. Running H1 K=300, 1000.

### 2026-10-05 destination-abundance bulk transport diagnostic

Fixed q_destination*(percell_target/m_destination)^(ab/ac=.4), normalized, mixed with original pooled profile.25% blend improves desired native-axis MSE all3 slightly; retrieval improves H1 but worsens Jurkat/HepG2. Raw ratio without ac/ab calibration over-amplifies and is not the tested calibrated recipe. No held-grid selection, generated cells, or full scorer. Not promoted. src/bulk_expression_transport.py; data/calibration/bulk_expression_transport.json.
  - **Result: topk loses badly.** H1 tk300 .178 / tk1000 .186 vs u3b15 .242. Calls did shrink (median 4,644 -> 160 / 342)
    but sign accuracy of calls stayed .55 and the kept genes are not the real DE ones (jac .093 -> .019). Raw fid
    .529 -> .308. With n_pred >> n_real, F = k/n_pred ~ sign accuracy over all calls; shrinking below n_real caps F at
    n_pred/n_real. Lesson: our most-confident genes are not the real DE genes on H1; set size is not the lever,
    which-genes is. Dropped.

### 2026-10-05 generated-cell variance audit

Eight fixed random H1 targets: active-gene generated/real logCP10k variance median 0.3846544952905322; real/control median 1.001725198846634. Seed1 variance with seed0 false-call groups; descriptive only. data/calibration/h1_generation_variance_audit.json. Current-combo pool1 full test queued after paired seed comparison, max2-active gate. Existing validated package stays pool4.

### 2026-10-05 paired seed result and pool1 variance

H1 seed1 combo .243283 vs u3b15 .241544, fresh shared-reference verification; plain MSE improves significantly under paired-target bootstrap. H1 seed0 edge also positive. HepG2 matched bundle instead favors u3b15 .293999 vs combo .270404, so no universal candidate winner. Pool1 combo scorer now running; direct eight-target variance audit pred/real median1.034 vs pool4 .385. Existing validated pool4 package unchanged and not uploaded.

### Codex: completed pool1 DE diagnostic (2026-10-05)

- 2026-10-05 (Codex): Completed exact-cache pool DE audit (150targets/context; own gene excluded). H1 pool4->1: median calls4189.5->1049 vsreal268, mean recall.533->.271, precision.130->.160, intersect-sign accuracy.656->.714. Hep calls34->16 vsreal54.5, recall.105->.072, precision.188->.252. Thus variance restoration trades away true response detection despite improving conditional precision/signs. Not official metric replacements. Evidence data/calibration/combo_pool_de_diagnostic.json. Avoid further global variance-only tuning; prioritize response-specific calibration learned on other contexts.


### Codex exact-DE presence ranking pilot

- 2026-10-05 (Codex): NEW exact-DE presence feasibility pilot completed, src/exact_de_presence_pilot.py. Trains on other two scored public contexts, target-disjoint validation; held labels reporting only, excludes own gene. Held average precision beats copied magnitude on all3: H1 .2711 vs .2230; Jurkat .1262 vs .0842; HepG2 .1473 vs .1266. AUC~.736/.749/.750. Absolute probability calibration fails (H1 predicts1.57% vs11.30% real, Jurkat7.67% vs3.67% real). Promising as within-target ranking ONLY; do not use raw probabilities to set response counts. No cells/scorer/upload. Results data/calibration/exact_de_presence_pilot.json. Next: context-robust rank weighting with training-only selection and discrimination checks.
- dsupp=0 (zero per-cell change on genes whose sources split in sign; median 10.8k/5.6k genes per target on H1/Jurkat),
  on u3b15: H1 .218 / Jurkat .221 vs .242 / .241. Dropped. Pattern from topk, dsupp, cons10: removing per-cell change
  always loses; boosting confident genes (u3b15/u3b2) wins. The per-cell lever is "add power where sure", not "de-call".
- u4b2 (cmin 4, x2): H1 .241 / Jurkat .218 vs u3b15 .242/.241. Boosting fewer genes loses; breadth matters. Next: cmin 2.
- u2b15 (cmin 2, x1.5; 6.7k/2.5k genes boosted): H1 .239 / Jurkat .237 vs u3b15 .242/.241. cmin 3 is the sweet spot. u3b15 stays first for 05:31 IST.
- **2026-10-06 u3b15 leaderboard .1920 (rank 238)** vs rmpc5 .1943. pds .716 (+.003), mse .168 (+.007), nmae .130
  (-.037), fid -.008 (+.003), reach .144 (+.010), jac .002. The boost buys reach/mse but costs fold-change size (nmae),
  exactly H1's pattern (tie, nmae -.032 -> -.051); Jurkat (+.031) overstated it. H1 again is the better board proxy.
  u3b2 (H1 nmae -.077) not uploaded. Idea: keep the boost's p-value effect but undo its fold-change inflation.
- u3n (u3b15 + cnorm: boost but keep per-target per-cell norm): H1 .241 (tie rmpc5), Jurkat .222. Reach gain vanishes with the inflation. Per-cell amplitude levers exhausted.
- Predicted-PC removal on rmpc5 (ppk=K ppg=-.5: remove half of the top-K principal directions of the panel's predicted pooled changes, renorm). Offline retrieval .858 -> .864 (K2) / .868 (K5). H1 full: pp2 .2421 (pds .778 -> .789, mse .073 -> .069), pp5 .241. Small (+.0013 H1, maybe +.002 board).
- Expression-dependent transfer (scratch expr_transfer.py): per source, sign accuracy on held-out real-DE pairs drops when dest expression >= 4x source (log2 ratio >= 2): Jurkat hct116 .66 -> .47, hek293t .65 -> .51; HepG2 hct116 .71 -> .53; H1 weak (.61 -> .56 k562). ~10% of pairs. Candidate per-source per-gene down-weight in fusion; expected gain small.
- xg1 (per-source per-gene weight min(1, 4(src+1)/(dest+1)) in fusion) on pp2: H1 .241 vs pp2 .242, Jurkat .210 vs rmpc5 .210 (mse down). Tie; dropped. Option kept (A.EXPR_W / xg).
- Layer-1 check: LINCS (Codex caches) as sign source on top-50 ranked real-DE pairs: LINCS sign acc .51-.60 (shRNA top-quartile |z| .55-.71), coverage 17-33% of pairs; ours when shRNA agrees .711/.802/.878 vs disagrees .671/.752/.801 (H1/Jurkat/HepG2). Consensus (not per-line) data, so no per-line transfer-rule training possible; weak voter only. pp2 not uploaded (user: gain too small).

## 2026-10-06: Layer 0, context identity (src/line_identity.py, DepMap 24Q4 from figshare 27993248, CC BY 4.0)
Control pseudobulk log2(CPM+1) vs DepMap log2(TPM+1), centred by the DepMap gene mean, top-2000 variable genes, Pearson.
Sanity: Jurkat -> JURKAT .85, HepG2 -> Hep G2 .80, K562 -> K-562 .67, HCT116 -> HCT 116 .87, RPE1 -> RPE1 .83.
- **ctx_A -> JURKAT .89** (next PF-382 .79): T-ALL, very likely Jurkat.
- **ctx_B -> HeLa .74** (next .40): almost certainly HeLa.
- **ctx_C -> CAL-33 .78** (next CAL 27 .67): head & neck squamous, likely CAL-33.
- The 300 2026 targets have **0** overlap with Nadig's 2,394 Jurkat/HepG2 targets (no direct copy).
- Consequence: our Jurkat local eval is context A's own line (different platform). Test ctxA_u3b2: rmpc5 on B/C,
  rmpc5 + u3b2 boost on A only (`A:cq=1 A:cboost=2 A:csupp=1 A:cmin=3`, new per-context override syntax).
- DepMap CRISPRGeneEffect: target essentiality predicts response size (Jurkat non-essential median 13 DE genes vs
  122-227 essential), but own-line effect is no better than K562 or the pan-line mean (spearman -.22/-.26/-.24 Jurkat).
  **2026 panel is ~all non-essential** (Jurkat: 14/300 < -.5, 0 < -1; CAL-33: 9/300, 0), unlike our local evals.
  Raw per-target metrics on non-essential targets only (scratch subset.py): H1 u3b15 raises nmae .961 -> .975 (matches
  the board), fid flat; Jurkat u3b2: fid .35 -> .46, reach .40 -> .46, nmae .867 -> .847, but jac .118 -> .030.
  Use the non-essential subset as the board-like proxy from now on.
- **Option 3 (learned transfer) feasibility, 2026-10-06** (scratch transfer_signal.py): data/lines has 9,526 targets
  common to K562 GW, HCT116 and HEK293T (X-Atlas). Per-target response correlation between lines (panel-mean removed,
  5,741 genes expressed in all 3, own gene out): non-essential (DepMap pan-mean > -.3, n=7,168) median r .021
  (HCT116~HEK293T) / .012 / .011 (with K562); moderately essential .03-.05; essential (< -1) .08-.10. The 2026 panel is
  ~all non-essential, so cross-line transfer carries almost no signal there: a learned transfer model has little to
  learn. Explains the plateau of all copy/vote/neural approaches.
- Within-line split-half (scratch splithalf.py), panel-mean removed, own+panel genes out: H1 non-essential targets
  median r .325 (76% > .2) at ~100 cells/half; Jurkat non-essential .17. So weak knockdowns have a REAL, reproducible,
  target-specific response, but it is line-specific (cross-line r ~.02 above).
- Control co-expression prior (pred = -corr(gene, target) over 5,000 control cells; scratch coexpr.py): correlation with
  the real target-specific change ~0 (H1 .004, Jurkat -.013) vs ours .055 / .112. Dead end.
- Net: ceiling for weak targets is large (split-half .33 vs ours .055 on H1) but no zero-shot source found that carries
  the line-specific part. Option 3 (learned transfer) has nothing to learn from for the 2026 regime.

## 2026-10-06: Bet 2, pretrained models
- X-Cell (Xaira): weights not released (HF repo README-only). Stack (Arc): in-context, needs perturbed prompt cells.
- Arc State, st-se-replogle-full/jurkat_0.99 (held-out Jurkat), published predictions (data/state/jurkat_pred.h5ad):
  head-to-head on 73 targets shared with our Jurkat eval (data/calibration/state_vs_ours.py):
  ours rmpc5 corr .125 / retrieval .928 / sign_top50 .827; State .154 / .829 / .821. Non-essential (15): ours
  .117/.927/.807, State .187/.891/.725. Comparable, not a big jump.
- **All public State Replogle checkpoints use one-hot perturbations over 2,024 essential-screen genes: 0 of the 300
  2026 targets are in the vocabulary.** Using State for 2026 needs retraining with gene-embedding perturbation
  features on genome-wide data (K562 GW + X-Atlas): a GPU project, still bounded by weak cross-line transfer.
- Realistic cell variance (pool=1/2 on rmpc5), H1: pool2 .238, pool1 .229 vs .241. Non-essential subset raw: pool1 plain mse .0020 -> .0015 (real gain), reach .18 -> .21, pds .883 -> .886, but fid .52 -> .46 and jac .073 -> .049. Raw fid there (~.52) is near the board's raw level, so the fid loss likely transfers; net expected negative. Not uploaded.
- **PIE (Arc, 2026-10-05; replogle_xdataset checkpoint, trained on Tahoe/Jiang/VCC25/X-Atlas, no Jurkat labels)** run on
  Modal (src/modal_pie.py; ~$1-2) for our 150 held-out Jurkat targets (data/calibration/pie_vs_ours.py):
  - alone, much weaker than ours: corr -.015 vs .060, retrieval .550 vs .844, sign_top50 .725 vs .781.
  - but complementary on direction: real-DE pairs, ours right .703; PIE agrees .801 / disagrees .555.
    Within >=3-source-unanimous genes: PIE agrees .902 (n=9,385) / disagrees .720; not unanimous: .752 / .516.
  - "unanimous & PIE agrees" is the first gene set at the 90% purity reach needs. Candidate: target the u3 boost
    (or p-value order) at that set only. Caveats: H1 can't test it (VCC25 H1 is in PIE's training); HeLa/CAL-33
    have no PIE context text (needs an OpenAI embedding or a stand-in context).
  - PIE-gated boost (u3 & PIE-agree only, ~944 genes/target; pie=1): Jurkat u3pie1.5 .232 vs u3b15 .241; u3pie2 .245 vs u3b2 .256. Narrower set loses (fid/reach down). PIE adds no usable lever here. Modal spend ~$3.
- **Why our board mse lags (cell_eval2 0.16 source + H1 breakdown, 2026-10-06):** the noise credit
  (r * min(tr Sigma_pred/n, k tr Sigma_real/n), k=1, #348) is only ~3% of our error (H1 capped .00299 vs
  uncapped .00290), so the gap is real accuracy. The score is sum(error)/sum(real change): on H1 the 20 strongest
  knockdowns are 67% of it and our ratio there is .916. For those 20, cosine(fused source change, real) median
  .20 (rest .03) and the LS-optimal scale on the raw fused change is .40: the pattern, not the size, is wrong,
  so amplitude tuning caps out (~cos^2). Several are stem-specific TFs (PRDM14, SOX2, SALL4) with no
  cross-line analogue. Top teams' mse .25-.70 imply much better patterns on the strong knockdowns.
- **State retrained (Modal, src/modal_state.py; ESM2 perturbation features; K562 GW + RPE1 + Nadig Jurkat/HepG2,
  <=40 cells/target, ~708k cells, 2,000 HVGs of 6,124 shared genes; 20k steps on L40S, ~45 min), H1 held out
  zero-shot:** chance level. On 150 H1 targets (data/calibration/state_pilot_vs_ours_h1.py): corr .004, retrieval
  .514, sign_top50 .509 vs ours rmpc5 .108 / .821 / .661. State's own eval: pearson_delta -.001, discrimination
  cosine .513. Predictions vary per target but in random directions (median inter-target cosine .055). Dropped.
  Modal spend for the pilot ~$6-8 (of $30 free credit).
