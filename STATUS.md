# Status (2026-09-30)

## Latest leaderboard (2026-10-03): best **0.1877** (aa75_c125, rank 201); board shows 0.1860 (aa75_c15). Goal: top 100 (≈0.215).

## Goal

Raise the leaderboard score above our current entry, 0.092 (α = 1). Our best so far was 0.097 (α = 0.5).

Rules we agreed on:
- Submit only real improvements. The leaderboard keeps only the newest entry, so a worse submission replaces a better one.
- A candidate must beat **0.178 on the local H1 test** first. H1 tracks the leaderboard; HepG2 doesn't.
- Ask before every `vcc submit`.

## Latest result

The "remove the shared part" idea (γ = 0) did **not** help on H1:
α = 0.5 → 0.173 and α = 1 → 0.176, against the 0.178 bar. Nothing was submitted.

## Where things stand (2026-10-01 early morning)

Nothing beat the 0.178 bar on H1, so nothing was submitted. Tried and dropped:

| idea | H1 |
|---|---|
| remove the part shared by all targets (γ = 0) | 0.173 / 0.176 |
| less random spread in the sampled cells (φ = 0) | 0.179 (tie; seed noise ~0.001) |
| **add H1's typical response** (MODEL.md) | 0.162 / 0.134 (targets look alike, pds drops) |
| average each K562 profile with its similar knockdowns | worse on the fast pds check |
| zero out weak genes / cap big changes | worse or tie on the fast pds check |

What I learned:
- The score is mostly pds (telling targets apart), and pds comes only from each target's own K562 profile.
  Anything added to every target hurts it.
- Every way of reshaping K562's profiles has been a tie or a loss. We're at the limit of what K562 alone gives.
- RPE1, HepG2 and Jurkat cover none of the 2026 targets; only K562 (272) and H1 (25) do.

## Latest (2026-10-01)

- Neighbor-gene rule: +0.002–0.003 on H1 in all 3 settings (0.180 vs 0.177/0.178). Kept in the model.
  The candidate `nb_a05` was built and validated but not submitted: the gain is too small (your call).
- X-Atlas (HCT116, HEK293T; 127 GB streamed): no gain. Both lines are much worse than K562 at predicting H1,
  and adding them to K562 is noise. Knockdown effects are mostly cell-type specific. Details in EXPERIMENTS.md.
- Conclusion: more unrelated cell types won't help. What would help is a knockdown screen in a cell type like the
  2026 contexts (stem-cell-like, like H1).

## Running now: the neighbor-gene effect

Switching a gene off also lowers the genes that start right next to it on the DNA. This holds in every cell type:

| genes starting near the target's start | H1 | K562 |
|---|---|---|
| within 1 kb | 0.28× (91% drop by >30%) | 0.60× |
| within 5 kb | 0.40× | 0.63× |
| random genes | 1.0× | 1.0× |

Rule (`add_neighbours` in `src/predict_2026.py`; gene positions from UCSC refGene, `data/annot/tss.csv`):

    Δ_{t,g} = min(Δ_{t,g}, d)      if g starts within 1 kb of t     (d = typical own drop ≈ log 0.31)
    Δ_{t,g} = min(Δ_{t,g}, d/2)    if g starts within 1–5 kb of t

- It affects 54 of the 300 2026 targets, including some of the 28 targets no line covers.
- On the fast pds check it helped a little (0.124 → 0.118 share of wrong targets ranked closer, lower is better).
  Those neighbors are also strongly changed genes in the truth, which should help the significant-gene scores.
- H1 runs: `nb_a1` (α = 1, γ = 0) and `nb_a1_g1` (α = 1). Bar: 0.178.

Other ideas (new information, not reshaping):
1. For the 25 targets H1 covers, use H1's profile alone (H1 looks like the contexts; K562 doesn't).
   That's a small number of targets, and it can't be tested locally.
2. Side effects on neighboring genes: switching off a gene often also lowers the gene next to it on the DNA.
   That holds in every cell type, including for the 28 targets no line covers. It needs gene positions,
   a small public annotation file.
3. Another large public screen covering more of the 300 targets in a cell type closer to H1.

## Current model (short)

Average each target's log share change over the public lines where that gene was knocked down. Scale it by α. Set the knocked-down gene to about 0.31× its normal level. Apply that to each context's control profile. Then sample 400 cells per target with negative-binomial noise fitted on that context's controls.
