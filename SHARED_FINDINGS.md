# Shared findings (Claude + Codex)

One place for what both agents have learned, so neither repeats the other's dead ends.
Keep entries short and add the date. Details live in EXPERIMENTS.md (Claude) and CALIBRATION.md,
SOURCE_SELECTION.md, LEARNED_SHRINKAGE.md (Codex).

## Facts about the task (2026-10-01)

- The 2026 contexts are **cancer cell lines from different tissues** (Arc, Cell 2026). Arc chose the targets
  "to provide a strong set of perturbations and responses". The contexts' controls look most like H1
  (corr 0.50–0.71 vs K562 0.22–0.33), probably because they share Arc's lab and pipeline, not their biology.
- **Coverage of the 300 targets:** K562 272, H1 25, RPE1/HepG2/Jurkat **0**, X-Atlas HCT116/HEK293T 300.
  So without X-Atlas, a 2026 prediction is K562 for 247 targets. Source weighting and HepG2 results don't
  reach those targets.
- **H1 is the local test that tracks the leaderboard;** HepG2 does not (alpha 1: HepG2 0.374 vs 0.171,
  leaderboard 0.092 vs 0.097, H1 0.178 vs 0.177).
- **Noise on the H1 test:** a different sampling seed moves the score by ~0.001. On HepG2 it's ~0.015.
- **Scorer:** pds and mse compare log1p(50000 · gene sum / total sum); pds excludes all panel target genes;
  DE uses a Wilcoxon test on per-cell normalized counts (Codex).
- The leaderboard keeps only the newest submission. The user approves every submit and only wants clearly
  better ones.

## Small real gains (H1, relative to baseline; bar alpha 0.5 = 0.177, alpha 1 = 0.178)

| change | H1 | who |
|---|---|---|
| vetted K562 control guides (514 core) | 0.1787 vs 0.1767 | Codex |
| neighboring genes drop too (start within 1 kb: own drop; 1–5 kb: half) | 0.180 vs 0.177 / 0.178, 3 of 3 settings | Claude |

These are independent, so they should stack. Neither is big enough to submit alone.

## Dead ends (don't repeat)

| idea | result | who |
|---|---|---|
| any shrinking below alpha 0.5; denoise per gene (z²/(z²+k)); per-target strength weighting | worse | both |
| alpha 2 | 0.069 on H1 | Claude |
| remove the panel-shared change (gamma 0) | tie | Claude |
| add H1's typical response to every target | 0.162 / 0.134; pds collapses | Claude |
| average a K562 profile with its similar knockdowns | targets harder to tell apart | Claude |
| zero out weak genes / cap big changes | tie or worse | Claude |
| K562 → H1 ridge map added to K562 | −0.008 [−0.026, +0.008] on the fast check, noise | Claude |
| learned per-target shrinkage (ridge) | extrapolates badly, H1 multipliers hit the bound | Codex |
| PCA + ridge learned transfer | hurts pds on H1/HepG2/Jurkat | Codex |
| depth-aware / factor noise in generated cells; plain Poisson | within noise | both |
| logbulk effect space; batch-corrected K562 | worse or tie | Codex |
| X-Atlas HCT116/HEK293T as extra sources | much worse than K562 alone at predicting H1, HCT116 or HEK293T; adding them to K562 is noise | Claude |
| batch-matched controls for X-Atlas | no change in split-half agreement | Claude |
| gradient boosting cross-line model (src/gbm.py) | H1 0.089 / 0.162 (rescaled) vs 0.180; closer on average but blander, pds and fid drop | Claude |

## In progress

- Claude: stacking vetted K562 controls (`core=1`, Codex's k562_core.npz) + neighbor rule on H1.
- Codex: similarity-weighted sources (HepG2 0.197 vs 0.188); H1 confirmation queued.
  Note that for 2026 this only affects the 25 H1-covered targets, unless X-Atlas lines are added as sources.

## Suggested next steps

1. Stack vetted K562 controls + the neighbor rule and test on H1 (two small independent gains).
2. If the GBM wins on H1, combine it with vetted controls.

## Notes between agents

Add a dated line under your own heading. Read the other heading before starting a new experiment.

### Claude → Codex

- 2026-10-01: Please read "Coverage of the 300 targets" above before tuning on HepG2. For 2026, 247 targets come
  from K562 alone, so HepG2-only gains and source weighting mostly don't reach the leaderboard.
  Judge candidates on H1.
- 2026-10-01: The H1 seed noise is ~0.001. A 0.002 gain on H1 is real but small; the user only wants
  submissions with clearly bigger gains.
- 2026-10-01: Your vetted-control change and my neighbor rule should stack. Could `calibration.py` expose the
  core-control K562 source so `src/predict_2026.py` / `src/local_eval.py` can use it (e.g. a `k562_core`
  entry in SOURCES)? I'll test the stack on H1.
- 2026-10-01: Shared data: `data/lines/{hct116,hek293t}.npz` (X-Atlas, all 300 targets) and
  `data/annot/tss.csv` (gene start positions) are free to use.
- 2026-10-01: We both run scorer jobs on this laptop (18 GB RAM, swap often full). Please keep to 2 jobs at once
  between us. Check `ps aux | grep cell-eval2` first.

- 2026-10-01 **IMPORTANT, split of work.** Public #82 entry (0.155 vs our 0.097) is open source (MIT):
  https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 (model.py, predict.py). Their recipe:
  1. Sources (weight): K562 per-cell CPM statistics (2), H1 (2), X-Atlas HCT116 (1), HEK293T (1) with batch-matched
     controls and a 1e5-count prior; CD4 T-cell genome-wide CRISPRi DE stats (Marson 2025, GWCD4i.DE_stats.h5ad, 0.5).
     Each source has its panel-shared change removed (common_subtract=1).
  2. **Two separate targets per knockdown:** the per-cell-normalized mean (log2 CPM ratio, amplitude 0.6) and the pooled
     profile (log1p(50000·share) delta, amplitude **0.3**), both clipped at ±3. Promoter-neighbor prior: genes
     within 500 bp keep 15%, ramping back to normal by 5 kb.
  3. **Cell generator `dual_moment_counts`:** 400 templates, each the average of 4 random real control cells, reweighted
     so the group's per-cell-normalized mean AND pooled profile hit the two targets exactly; then integer rounding
     that keeps each cell's total. No gamma-Poisson sampling noise.
  **Claude does:** (2) + (3) on our current sources, tested on H1 (`src/atlas_shift.py`, H1 runs named `as_*`).
  **Suggested for Codex (parallel, no overlap):** the source side, in our `data/lines/{name}.npz` format so both
  pipelines can use it:
  (a) CD4 T-cell source from GWCD4i.DE_stats.h5ad (16.8 GB, ask the user before downloading);
  (b) K562 per-cell mean-CPM statistics (needs K562_gwps_raw_singlecell, 66 GB; ask the user), or check whether our
      bulk-share approximation is close enough;
  (c) X-Atlas with batch-matched controls + count prior (my `src/xatlas.py` uses global controls; batch-matching
      didn't change split-half agreement in my 12-batch test, but #82 uses it).
  Please note in "Codex → Claude" which of these you take, so we don't both do it.

- 2026-10-01: We both launched "vetted controls + neighbor rule, alpha 0.5" on H1 (yours: codex_core_neighbours_a05).
  I stopped my copy and kept alpha 1 (`core_nb_a1`), so yours covers alpha 0.5. Please check this file before launching.
- 2026-10-01: #82 generator running on H1 as `as_c06_b03` (src/atlas_shift.py, vendored their model.py under
  src/third_party/, MIT). Dry run: hits both targets within 0.2%, 0.3 s per target.

- 2026-10-01 **RESULT:** #82 recipe on our sources scores **0.205 on H1** (vs 0.180 best before; +0.025 ≈ 25x seed noise).
  pds 0.71, mse 0.009 (first time above 0), nmae −0.02, fid 0.44. Vetted controls + neighbor rule at alpha 1: 0.178 (no stack).
  The generator + two separate scales are the win, so tune and build on `src/atlas_shift.py` from here on, not predict_2026.py.
  Building the 2026 file `as_c06_b03` now (user approves submits). Testing ab=0.5 on H1 (`as_c06_b05`).
  Most useful Codex work now: the missing sources (CD4 T cells; K562 per-cell CPM; batch-matched X-Atlas) as
  `data/atlas_shift/{line}.npz` (keys targets, ec, eb: centered changes on the 2026 gene list, see `build()`).
  Then add the line to `WEIGHTS` in atlas_shift.py and test with `python src/atlas_shift.py local h1 NAME`.

### Codex → Claude

(empty)
