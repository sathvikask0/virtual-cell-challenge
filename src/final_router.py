"""Pick a recipe per 2026 context from its control cells (for the final D/E/F, and as a check on A/B/C).

Arc's dataset note says the 2026 lines differ from the 2025 H1 ESCs. Its cell-type illustration spans both years;
it does not establish H1 or another stem context among D/E/F. This router is an experimental control-based rule,
not a prediction of final-set identities. Stem-source gains on H1 do not establish transfer to unseen final lines.

Rule: stem score = mean log1p(CPM) of pluripotency markers in the context's controls. Above STEM_MIN -> "stem" recipe.

Usage (from src/): python final_router.py            (prints each context's markers and chosen recipe)
"""
import numpy as np

import atlas_shift as A

MARKERS = ["POU5F1", "NANOG", "LIN28A", "DPPA4", "TDGF1", "L1TD1"]
STEM_MIN = 3.0  # mean log1p(CPM); H1 controls ~6-7, A/B/C ~0-0.5
RECIPES = {
    "default": dict(alpha=.75, beta=.75, ac=1.25, cd4=1.0, pcs=100, pcg=-1, pchi=5, pcn=1),  # rmpc5, leaderboard .1943
    "stem": dict(alpha=.75, beta=.75, ac=1.25, cd4=1.0, kolf=2.0, ipsc=2.0),  # H1 .2712 vs .2365 default; EXPERIMENTS.md
}


def stem_score(genes, ctrl):
    cpm = 1e6 * np.asarray(ctrl, float) / np.sum(ctrl)
    gi = {g: i for i, g in enumerate(genes)}
    v = [np.log1p(cpm[gi[g]]) for g in MARKERS if g in gi]
    return float(np.mean(v)) if v else 0.0


def route(genes, ctrl):
    s = stem_score(genes, ctrl)
    return ("stem" if s >= STEM_MIN else "default"), s


if __name__ == "__main__":
    for name in ["ctx_A", "ctx_B", "ctx_C", "h1", "k562", "jurkat", "hepg2"]:
        z = np.load(A.ROOT / "data/lines" / f"{name}.npz", allow_pickle=True)
        r, s = route(z["genes"].astype(str), z["ctrl"])
        print(f"{name:7s} stem score {s:5.2f} -> {r}")
