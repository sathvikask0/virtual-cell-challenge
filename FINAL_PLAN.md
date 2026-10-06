# Final round plan (contexts D/E/F)

Final test set released **Thu Oct 22**; submissions due **Thu Nov 5, 23:59 UTC**. Arc: "Participants may use any
modeling strategy and train their models on any data." Open questions to check on Oct 22 (site is JS-only, check in a
browser or with `vcc`): number of final submissions, which one counts (last / best / selected), whether final scores
are shown before the deadline.

## Day-of steps (one command after the download)

1. Download (needs the user's OK): `.venv/bin/vcc datasets list`, then `.venv/bin/vcc datasets download <id>`; unzip to
   `data/vcc/final/` (expect `context_D.h5ad`, `context_E.h5ad`, `context_F.h5ad`, `gene_names.csv`,
   `pert_counts.csv`).
2. Run the pipeline with the current best recipe (rmpc5):

   ```
   ./final_pipeline.sh data/vcc/final "D,E,F" final_rmpc5 alpha=.75 beta=.75 ac=1.25 pcs=100 pcg=-1 pchi=5 pcn=1
   ```

   It builds each context's control profile, prints the DepMap identity of each context, rebuilds every source cache
   for the final target panel into `data/atlas_shift_final_rmpc5/` (old caches untouched), then builds and validates
   `data/submissions/final_rmpc5.vcc`. It never uploads.
3. Read the identity printout and apply the recipe table below (a second build with per-context options).
4. Upload only with the user's explicit OK.

Rehearsal on A/B/C (`dry_rmpc5`, 2026-10-06): all steps ran in **12 minutes** and the output was **byte-identical** to the
uploaded `rmpc5.vcc` (log: `data/calibration/dry_pipeline.log`). Upload takes ~20-30 more minutes.

## Recipe by identity (from src/line_identity.py; trust a match at r >= 0.8 with a clear gap to the next line)

| Identity of a final context | Recipe | Why |
|---|---|---|
| **HCT116, HEK293T or K562** (we have that line's own genome-wide screen) | rmpc5 with that source weighted very high, e.g. `hct116=20` | Within-line responses repeat far better than cross-line (split-half r ~.33 vs ~.02 on weak knockdowns). Potentially the biggest single gain. Use `D:hct116=20` (per-context, verified). |
| **Jurkat / T-cell** | rmpc5 (u3b2 on that context only is the Jurkat-local winner but did not hold on the board as u3b15) | ctxA_u3b2 never tested on the board. |
| **Stem-like** (src/final_router.py stem score >= 3) | stem recipe: + `kolf=2 ipsc=2` | H1 .2712 vs .2365. |
| anything else | rmpc5 | Board best .1943. |

## TODO before Oct 22
- Per-context source weights (`D:hct116=20` etc. in `agree_alloc.py`) verified 2026-10-06: a full build with `A:k562=10` validated and changed only context A.
- Measure the same-line gain locally: e.g. HCT116 as destination is not possible (no local real HCT116 cells), so use
  Jurkat/HepG2 split halves: predict half B from half A of the same line vs from other lines.
