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

---

## Current

### Exp 4: practice on unseen cell lines (leave one cell line out)
- Goal: copy the 2026 setup. Fit on some cell lines, predict another from its control cells only.
- Data (public, free):
  - K562 genome-wide screen (Replogle 2022, averaged file): 9,866 targets
  - RPE1 essential screen (Replogle 2022, averaged file): 2,393 targets
  - HepG2 and Jurkat essential screens (Nadig 2025, per-cell files, ~15 GB): downloading
  - H1 (VCC 2025): 300 targets
- `src/lines.py` puts every line in one format; `src/cross_line.py` runs the test.
- Scored on the 6,633 genes measured in every line (the Replogle files keep only ~8k genes).
- 272 of the 300 2026 targets were switched off in K562. None were in the essential screens.
- Predictors: control (no change), mean (average change), mean+drop (switched-off gene itself drops),
  transfer (reuse the same gene's change from another line).

First run (H1, K562, RPE1), same metrics as before; "seen" = gene switched off in a training line:

| held out | targets | predictor | pearson_delta | mae | discrimination |
|---|---|---|---|---|---|
| h1 | seen (280) | control | 0.00 | 0.402 | 0.500 |
| h1 | seen (280) | mean+drop | 0.09 | 0.606 | 0.492 |
| h1 | seen (280) | **transfer** | **0.21** | 0.745 | **0.353** |
| k562 | seen (300) | mean+drop | −0.18 | 0.145 | 0.498 |
| k562 | seen (300) | **transfer** | **0.15** | 0.251 | **0.335** |
| rpe1 | seen (300) | mean+drop | 0.08 | 0.277 | 0.499 |
| rpe1 | seen (300) | **transfer** | **0.16** | 0.300 | **0.435** |
| k562 | unseen (300) | mean+drop | 0.01 | 0.106 | 0.498 |

- So far: reusing the same gene's change from another line clearly helps correlation and discrimination.
- But it makes mae worse, and the average change is worse than "no change" on mae in every line.
  Changes copied across lines carry a lot of noise.
