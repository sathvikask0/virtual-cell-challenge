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

