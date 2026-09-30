# Status (2026-09-30)

## Goal

Raise the leaderboard score above our current entry, 0.092 (α = 1). Our best so far was 0.097 (α = 0.5).

Rules we agreed on:
- Submit only real improvements. The leaderboard keeps only the newest entry, so a worse submission replaces a better one.
- A candidate must beat **0.178 on the local H1 test** first. H1 tracks the leaderboard; HepG2 doesn't.
- Ask before every `vcc submit`.

## Latest result

The "remove the shared part" idea (γ = 0) did **not** help on H1:
α = 0.5 → 0.173 and α = 1 → 0.176, against the 0.178 bar. Nothing was submitted.

## Running now (H1, ~45–60 min)

- `hc1_a1` and `hc1_a05`: the new model (H1 typical response + K562 specific part), at α = 1 and α = 0.5.

Done: φ = 0 scored 0.179 and a different seed scored 0.177, against 0.178. On H1 the seed moves the score by only ~0.001, so all our older ideas were ties.

## Most promising idea: use H1's typical response (next run: `hc`)

New finding:
- RPE1, HepG2 and Jurkat cover **none** of the 300 2026 targets. 247 are covered by K562 only,
  25 by K562 + H1, and 28 by no line.
- The 2026 contexts' control cells look much more like H1 than K562:

| | H1 | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|---|
| context A | **0.50** | 0.33 | 0.29 | 0.29 | 0.44 |
| context B | **0.63** | 0.22 | 0.29 | 0.25 | 0.22 |
| context C | **0.71** | 0.31 | 0.36 | 0.37 | 0.29 |

(correlation of log expression, 6,123 shared genes)

So today the "typical knockdown response" in our predictions comes from K562, the wrong cell type.

New model, per target t:

    Δ_t = μ_H1 + α·(Δ_K562,t − μ_K562),   own gene's drop kept

- μ = a line's average change over all its knockdowns (its typical response).
- Δ_K562,t − μ_K562 = what's specific to target t, taken from K562.
- The total-count change is also taken from H1's typical response.

Local test (fair): μ_H1 comes only from the 150 H1 knockdowns that are **not** in the 150 test targets.
It correlates 0.91 with the test targets' typical response; K562's correlates only 0.22.
Caveat: the local test is optimistic, because H1 → H1 is easier than H1 → 2026 contexts.

Command: `python src/local_eval.py predict h1 NAME alpha=1 gamma=0 hc=1`
Code: `line_common` and `add_common` in `src/predict_2026.py`.

## Current model (short)

Average each target's log share change over the public lines where that gene was knocked down. Scale it by α. Set the knocked-down gene to about 0.31× its normal level. Apply that to each context's control profile. Then sample 400 cells per target with negative-binomial noise fitted on that context's controls.
