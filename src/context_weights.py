"""Do source weights depend on the cell type being predicted? Test on Jurkat (T cell) and HepG2 (liver).

Why: the 2026 contexts are not stem cells (A is T-cell-like, B mesenchymal, C squamous epithelial), so source
weights tuned on H1 may be wrong. If a T-cell line wants much more CD4 than a liver line does, per-context
weights are worth it.

The atlas_shift / CD4 caches only hold the 2026 + H1 targets, so this rebuilds them, plus the Jurkat and HepG2
local-eval targets, into data/atlas_shift_x/ (the main caches are untouched).

Usage (from src/):
  python context_weights.py build
  python context_weights.py local LINE NAME k562=2 h1=2 hct116=1 hek293t=1 cd4=1 [ac=1 ab=.5]
  python context_weights.py submission NAME A=1 B=0 C=0 [ac=1 ab=.5]
  python context_weights.py submission NAME A=1 B=1 C=k562:2,h1:2,hct116:2,hek293t:1,cd4:0
                                            2026 build with a CD4 weight per context (main caches, default
                                            source weights); writes data/submissions/NAME.vcc, no upload
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np

import atlas_shift as A
import atlas_cd4 as C

XOUT = A.ROOT / "data/atlas_shift_x"
LINES = ("jurkat", "hepg2")


def local_targets():
    t = set()
    for line in LINES:
        p = A.ROOT / "data/local_eval" / line / "real.h5ad"
        if p.exists():
            t |= set(ad.read_h5ad(p, backed="r").obs["target"].astype(str))
    return t


def use_x():
    A.OUT = XOUT
    C.PATH = XOUT / "cd4_de.npz"
    base = A.wanted_targets
    A.wanted_targets = lambda: base() | local_targets()


def build():
    XOUT.mkdir(exist_ok=True)
    use_x()
    A.build()
    C.build()


def local(line, name, ac=1.0, ab=0.5, **w):
    use_x()
    cd4 = w.pop("cd4", 0.0)
    A.WEIGHTS = {k: float(v) for k, v in w.items() if float(v) > 0}
    print(f"{line}/{name}: weights {A.WEIGHTS}, cd4 {cd4}", flush=True)
    if cd4 > 0:
        A.profiles = C.make_profiles(cd4)
    A.local(line, name, ac=ac, ab=ab)


def with_weights(weights, cd4):
    """profiles function that uses these source weights (and CD4 weight) for one context."""
    base = C.make_profiles(cd4) if cd4 > 0 else A.profiles

    def profiles(*args, **kw):
        saved = A.WEIGHTS
        A.WEIGHTS = weights
        try:
            return base(*args, **kw)
        finally:
            A.WEIGHTS = saved
    return profiles


def submission(name, ac=1.0, ab=0.5, **spec):
    """spec: context -> CD4 weight (A=1), or context -> 'k562:2,h1:2,hct116:2,hek293t:1,cd4:0'."""
    for c, v in spec.items():
        if isinstance(v, str) and ":" in v:
            w = {k: float(x) for k, x in (kv.split(":") for kv in v.split(","))}
            cd4 = w.pop("cd4", 0.0)
            A.CTX_PROFILES[c] = with_weights({k: x for k, x in w.items() if x > 0}, cd4)
        else:
            A.CTX_PROFILES[c] = C.make_profiles(float(v)) if float(v) > 0 else A.profiles
    print(f"{name}: per context {spec}, default sources {A.WEIGHTS}", flush=True)
    A.build_2026(name, ac=ac, ab=ab)


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build()
    elif sys.argv[1] == "submission":
        kv = dict(a.split("=", 1) for a in sys.argv[3:])
        submission(sys.argv[2], ac=float(kv.pop("ac", 1.0)), ab=float(kv.pop("ab", 0.5)), **kv)
    else:
        opts = {k: float(v) for k, v in (a.split("=") for a in sys.argv[4:])}
        local(sys.argv[2], sys.argv[3], **opts)
