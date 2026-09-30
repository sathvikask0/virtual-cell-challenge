# Model math

Notation: g = gene, t = knocked-down target, ℓ = a public line (H1, K562, RPE1, HepG2, Jurkat),
c = a 2026 context (A/B/C), ε = 1e-5.

## 1. Change in one public line

For each line ℓ: x_{ℓ,t,g} = mean counts per cell for target t's cells, and x_{ℓ,0,g} = the same for control cells.
Shares: p = x / Σ_g x.

    Δ_{ℓ,t,g} = log( (p_{ℓ,t,g} + ε) / (p_{ℓ,0,g} + ε) )        change in gene g's share
    r_{ℓ,t}   = log( Σ_g x_{ℓ,t,g} / Σ_g x_{ℓ,0,g} )              change in total counts

## 2. A line's typical response

Average over the targets T_ℓ knocked down in line ℓ:

    μ_{ℓ,g} = (1/|T_ℓ|) Σ_{t∈T_ℓ} Δ_{ℓ,t,g}
    ρ_ℓ     = mean_t r_{ℓ,t}

## 3. Change predicted for target t

Δ̄_{t,g} = the average of Δ_{ℓ,t,g} over the lines that knocked t down.
For the 2026 targets that's K562 only for 247 targets, K562 + H1 for 25, and none for 28.
Targets knocked down in no line get the average change over all targets.

**Old model** (on the leaderboard: 0.097 at α = 0.5, 0.092 at α = 1):

    Δ_{t,g} = α · Δ̄_{t,g}
    r_t     = α · mean_ℓ r_{ℓ,t}

**New model** (being tested), where μ̄_g is the same average of μ_{ℓ,g} over the lines that knocked t down:

    Δ_{t,g} = μ_{H1,g} + α · (Δ̄_{t,g} − μ̄_g)
               └ typical ┘   └ specific to t ┘
    r_t     = ρ_{H1}

- The typical part comes from H1. The 2026 contexts' control cells correlate 0.50–0.71 with H1
  but only 0.22–0.33 with K562 (log expression, 6,123 shared genes).
- The specific part comes from K562, the only line that covers most 2026 targets.

In both models the target gene itself is set to d = median over ℓ,t of Δ_{ℓ,t,t} ≈ log 0.31 (not scaled).

## 4. Apply to context c

q_{c,g} = the context's control share.

    s_{t,g} ∝ max( (q_{c,g} + ε) · exp(Δ_{t,g}) − ε, 0 ),   normalised so Σ_g s_{t,g} = 1

## 5. Sample 400 cells per target and context

For each cell i:

    N_i ~ Uniform{ control totals of c },   N_i ← N_i · exp(r_t)
    w_{ig} ~ Gamma(shape 1/φ_g, scale φ_g)       (mean 1, variance φ_g)
    X_{ig} ~ Poisson( N_i · s_{t,g} · w_{ig} )

So each count is negative binomial with Var = μ + φ_g μ².
φ_g is fitted on the context's controls from the mean and variance of depth-normalised counts,
clipped to [0, 10].
Generating control cells this way and comparing them with real controls gives 4 false "changed" genes (real vs real: 0).

## Local test (H1 plays an unseen context)

- 150 random H1 targets, up to 200 cells each. H1's control cells are split in half: one half is the model's
  input, the other is scored with the answers.
- Scored with Arc's scorer (cell-eval2, preset vcc2026), relative to the official baseline.
- Δ̄ uses the other lines only (H1 is held out).
- μ_{H1} comes only from the 150 H1 targets **not** in the test set, which mimics 2026, where most targets are not in H1.
  It correlates 0.91 with the test targets' true typical response; K562's correlates 0.22.
  Caveat: this test is optimistic. Here the typical response goes H1 → H1; for real it goes H1 → the 2026 contexts.
- Bar to beat before submitting: 0.178 (old model, α = 1).

## Why it should help

The scorer's baseline is the true average response, so the per-gene error (mse) is already at its floor.
The significant-gene metrics (fid, reach, jac) depend mostly on the shared response plus the target's own drop.
We sit near or below the baseline on those because our shared part comes from K562.

## Code

- `src/predict_2026.py`: `source_changes` (Δ̄, μ̄), `target_change` (α, own drop), `recenter` (removes μ̄, γ = 0),
  `line_common` (μ_{H1}, ρ_{H1}), `add_common`, `context_stats` (q, totals, φ), `sample_cells`.
- Local test: `python src/local_eval.py predict h1 NAME alpha=1 gamma=0 hc=1`.
