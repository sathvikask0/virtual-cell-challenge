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

Next ideas (new information, not reshaping):
1. For the 25 targets H1 covers, use H1's profile alone (H1 looks like the contexts; K562 doesn't).
   That's a small number of targets, and it can't be tested locally.
2. Side effects on neighboring genes: switching off a gene often also lowers the gene next to it on the DNA.
   That holds in every cell type, including for the 28 targets no line covers. It needs gene positions,
   a small public annotation file.
3. Another large public screen covering more of the 300 targets in a cell type closer to H1.

## Current model (short)

Average each target's log share change over the public lines where that gene was knocked down. Scale it by α. Set the knocked-down gene to about 0.31× its normal level. Apply that to each context's control profile. Then sample 400 cells per target with negative-binomial noise fitted on that context's controls.
