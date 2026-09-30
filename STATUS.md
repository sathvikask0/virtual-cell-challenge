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

## What happens next

1. Try ideas that change *which* genes we predict to move, not just how much. The significance-overlap and direction scores (jac, fid) are near 0, so that's where the room is:
   - weight each public line by how clean its profile for that target is (few real changes = mostly noise)
   - predict changes only for genes with a clear signal in at least one line; set the rest to 0 but don't shrink the clear ones
2. Log results in `EXPERIMENTS.md` and push to `exp/calibration`.

## Current model (short)

Average each target's log share change over the public lines where that gene was knocked down. Scale it by α. Set the knocked-down gene to about 0.31× its normal level. Apply that to each context's control profile. Then sample 400 cells per target with negative-binomial noise fitted on that context's controls.
