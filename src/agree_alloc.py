"""Cross-source agreement allocation of each target's effect size (idea from leaderboard #22's description).

Per target t, with sources s (K562, H1, HCT116, HEK293T, CD4; per-cell-space changes E_c):
  A_t = mean pairwise cosine between the sources' centered changes (genes both measured)
  m_t = norm of the fused change
  s_t = max(A_t, floor)^alpha * m_t^-beta,   rescaled so sum_t s_t^2 |E_t|^2 = sum_t |E_t|^2   (energy kept),
        capped at smax (default 2)
Then the fused change of target t is multiplied by s_t: in the pooled target (bulk), and optionally the per-cell one.
Targets the sources agree on get bigger changes, conflicting ones shrink toward the control.
A target with fewer than 2 sources gets the median A.

Usage (from src/):
  python agree_alloc.py local LINE NAME [alpha=.75] [beta=.75] [cell=0] [cd4=1] [ac=1] [ab=.5] [tpow=1] [x=0]
  (tpow: Codex's template variance power; x=1: caches incl. Jurkat/HepG2 targets)
  python agree_alloc.py submission NAME [same options]
"""
import os
import sys

import numpy as np

import atlas_shift as A
import atlas_cd4 as C


def source_rows(targets, genes, exclude):
    """list of (targets x genes) E_c arrays, NaN where missing, one per source incl. CD4."""
    g26 = {g: i for i, g in enumerate(A.genes26())}
    cols = np.array([g26.get(g, -1) for g in genes])
    out = []
    for line in A.WEIGHTS:
        if line in exclude:
            continue
        d = np.load(A.OUT / f"{line}.npz")
        rows = {t: i for i, t in enumerate(d["targets"])}
        V = np.full((len(targets), len(genes)), np.nan, np.float32)
        for i, t in enumerate(targets):
            if t in rows:
                V[i] = np.where(cols >= 0, d["ec"][rows[t]][np.maximum(cols, 0)], np.nan)
        out.append(V)
    with np.load(C.PATH) as d:
        ci = {g: i for i, g in enumerate(d["genes"])}
        cc = np.array([ci.get(g, -1) for g in genes])
        rows = {t: i for i, t in enumerate(d["targets"])}
        V = np.full((len(targets), len(genes)), np.nan, np.float32)
        for i, t in enumerate(targets):
            if t in rows:
                V[i] = np.where(cc >= 0, d["ec"][rows[t]][np.maximum(cc, 0)], np.nan)
        out.append(V)
    return out


def agreement(rows, own):
    T = rows[0].shape[0]
    A_t = np.full(T, np.nan)
    for i in range(T):
        cos = []
        for a in range(len(rows)):
            for b in range(a + 1, len(rows)):
                x, y = rows[a][i], rows[b][i]
                m = np.isfinite(x) & np.isfinite(y) & ~own[i]
                if m.sum() > 100:
                    cos.append(x[m] @ y[m] / (np.linalg.norm(x[m]) * np.linalg.norm(y[m]) + 1e-12))
        if cos:
            A_t[i] = np.mean(cos)
    return np.where(np.isfinite(A_t), A_t, np.nanmedian(A_t))


def soft_generator():
    """Wrap the cell generator: if a target's two moments can't be fit, shrink both toward the templates' own
    moments by 10% steps until they can (counted in soft_generator.fallbacks)."""
    fit = A.dual_moment_counts

    def gen(template, probability, bulk_probability, **kw):
        t = np.asarray(template, float)
        base_c = (t / np.maximum(t.sum(1, keepdims=True), 1e-12)).mean(0)
        base_b = t.sum(0) / t.sum()
        for lam in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.0):
            try:
                out = fit(template, base_c + lam * (probability - base_c), base_b + lam * (bulk_probability - base_b), **kw)
                if lam < 1:
                    gen.fallbacks.append(lam)
                return out
            except ValueError as e:
                if "Moment fitting failed" not in str(e) and "projection too large" not in str(e):
                    raise
        raise RuntimeError("generator failed even with no change")
    gen.fallbacks = []
    return gen


def weighted_fusion(targets, genes, exclude, q, cd4, gamma, floor=0.002):
    """Fused (E_c, E_b) with per-target source weights W_s * max(c_s, floor)^gamma, where c_s is the cosine of
    source s with the mean of the other sources for that target (leave-one-out agreement)."""
    g26 = {g: i for i, g in enumerate(A.genes26())}
    cols = np.array([g26.get(g, -1) for g in genes])
    EC, EB, W = [], [], []
    for line, w in A.WEIGHTS.items():
        if line in exclude:
            continue
        d = np.load(A.OUT / f"{line}.npz")
        rows = {t: i for i, t in enumerate(d["targets"])}
        c = np.full((len(targets), len(genes)), np.nan, np.float32)
        b = c.copy()
        for i, t in enumerate(targets):
            if t in rows:
                c[i] = np.where(cols >= 0, d["ec"][rows[t]][np.maximum(cols, 0)], np.nan)
                b[i] = np.where(cols >= 0, d["eb"][rows[t]][np.maximum(cols, 0)], np.nan)
        EC.append(c); EB.append(b); W.append(w)
    if cd4 > 0:
        c = source_rows(targets, genes, exclude)[-1]
        base = 50000 * np.asarray(q)
        EC.append(c)
        EB.append(np.log1p(base[None, :] * np.exp2(np.clip(c, -10, 10))) - np.log1p(base[None, :]))
        W.append(cd4)
    own = np.array([genes == t for t in targets])
    S = len(EC)
    unit = [np.where(np.isfinite(x) & ~own, x, 0) for x in EC]
    unit = [u / (np.linalg.norm(u, axis=1, keepdims=True) + 1e-12) for u in unit]
    present = [np.isfinite(x).any(1) for x in EC]
    wt = np.zeros((S, len(targets)))
    for k in range(S):
        others = sum(unit[j] for j in range(S) if j != k)
        cs = (unit[k] * others).sum(1) / (np.linalg.norm(others, axis=1) + 1e-12)
        nother = sum(present[j] for j in range(S) if j != k)
        cs = np.where(nother > 0, cs, 1.0)  # a lone source keeps its weight
        wt[k] = np.where(present[k], W[k] * np.maximum(cs, floor) ** gamma, 0)
    num_c = sum(wt[k][:, None] * np.nan_to_num(EC[k]) for k in range(S))
    num_b = sum(wt[k][:, None] * np.nan_to_num(EB[k]) for k in range(S))
    den = sum(wt[k][:, None] * np.isfinite(EC[k]) for k in range(S))
    ec = np.where(den > 0, num_c / np.maximum(den, 1e-12), 0).astype(np.float32)
    eb = np.where(den > 0, num_b / np.maximum(den, 1e-12), 0).astype(np.float32)
    rel = wt / np.maximum(wt.sum(0, keepdims=True), 1e-12)
    print(f"  source weighting gamma={gamma}: mean weight share per source {np.round(rel.mean(1), 3)}", flush=True)
    return ec, eb


CTX = {}  # the current context's PCA basis of control-cell log expression, set by capture_controls()


def capture_controls(k=100, n=4000):
    """Wrap atlas_shift.control_stats so each context's control cells also give a PCA basis (genes x k) of
    log1p(10k-normalised expression), centred, from up to n cells."""
    from sklearn.decomposition import PCA
    import scipy.sparse as sp
    original = A.control_stats

    def stats(X, *a, **kw):
        out = original(X, *a, **kw)
        R = sp.csr_matrix(X, dtype=np.float64)
        idx = np.random.default_rng(0).choice(R.shape[0], min(n, R.shape[0]), replace=False)
        R = R[idx]
        lib = np.asarray(R.sum(1)).ravel()
        L = np.log1p((sp.diags(1e4 / lib) @ R).toarray())
        CTX["V"] = PCA(n_components=k, random_state=0).fit(L).components_.T  # genes x k
        print(f"  context PCA: {k} components from {len(idx)} control cells", flush=True)
        return out
    A.control_stats = stats


def pca_smooth(E, own, lam, renorm, g=0.0, lo=0, spec=0, hi=0):
    """lam: blend toward the projection, or g > 0: E + g * P. lo: skip the first lo PCs.
    spec: project only the target-specific part (E minus the mean over this panel's targets)."""
    V = CTX["V"][:, int(lo):(int(hi) or None)]
    E0 = np.where(own, 0, E)
    P = ((E0 - E0.mean(0) if spec else E0) @ V) @ V.T
    out = E + g * P if g else (1 - lam) * E + lam * P
    if renorm:
        n0 = np.linalg.norm(np.where(own, 0, E), axis=1, keepdims=True)
        n1 = np.linalg.norm(np.where(own, 0, out), axis=1, keepdims=True)
        out = out * n0 / np.maximum(n1, 1e-12)
    return np.where(own, E, out)


def make_profiles(alpha=0.75, beta=0.75, cell=0, cd4=1.0, floor=0.02, smax=2.0, srcw=0.0, nr=0.0, nrmax=3.0, gconf=0.0, ew=0.0, eb_boost=0.0, ek=500.0, pc=0.0, cq=0.0, cboost=1.6, csupp=0.5, ipsc=0.0, pcl=0.0, pcc=0.0, pcn=0.0, pcg=0.0, pclo=0, pcspec=0, pchi=0, cmin=2):
    def profiles(targets, genes, m, q, ac, ab, exclude=(), agree=0, thr=0.0):
        ec, eb = A.fused(targets, genes, exclude)
        den = np.zeros_like(ec)
        g26 = {g: i for i, g in enumerate(A.genes26())}
        cols = np.array([g26.get(g, -1) for g in genes])
        for line, w in A.WEIGHTS.items():
            if line in exclude:
                continue
            with np.load(A.OUT / f"{line}.npz") as d:
                ti = {t: i for i, t in enumerate(d["targets"])}
                for i, t in enumerate(targets):
                    if t in ti:
                        den[i] += w * ((cols >= 0) & np.isfinite(d["ec"][ti[t]][np.maximum(cols, 0)]))
        rows = source_rows(targets, genes, exclude)
        if srcw > 0:
            ec, eb = weighted_fusion(targets, genes, exclude, q, cd4, srcw)
        elif cd4 > 0:
            ec, eb = C.add_cd4(ec, eb, den, rows[-1], q, cd4)
        if ipsc > 0:  # Codex's iPSC 34-line LFC source (own-suppression QC), added like CD4 via destination controls
            with np.load(A.OUT / "ipsc_de.npz") as d:
                gi = {g: i for i, g in enumerate(d["genes"])}
                cc = np.array([gi.get(g, -1) for g in genes])
                ok = d["own_suppression_pass"]
                ti = {t: i for i, t in enumerate(d["targets"]) if ok[i]}
                V = np.full((len(targets), len(genes)), np.nan, np.float32)
                for i, t in enumerate(targets):
                    if t in ti:
                        V[i] = np.where(cc >= 0, d["ec"][ti[t]][np.maximum(cc, 0)], np.nan)
            den2 = den + (cd4 * np.isfinite(rows[-1]) if cd4 > 0 else 0)
            ec, eb = C.add_cd4(ec, eb, den2, V, q, ipsc)
            rows = rows + [V]
            print(f"  iPSC source: {np.isfinite(V).any(1).sum()} of {len(targets)} targets", flush=True)
        own = np.array([genes == t for t in targets])
        A_t = agreement(rows if cd4 > 0 else rows[:-1], own)
        E = np.where(own, 0, ec)
        m_t = np.linalg.norm(E, axis=1) + 1e-9
        s = np.maximum(A_t, floor) ** alpha * m_t ** -beta
        energy = (m_t ** 2).sum()
        for _ in range(20):  # energy kept, boost capped at smax (the cell generator can't fit much larger changes)
            s = np.minimum(s * np.sqrt(energy / ((s * m_t) ** 2).sum()), smax)
        profiles.scale = s
        profiles.A = A_t
        print(f"  agreement A: median {np.median(A_t):.3f}; scale s: min {s.min():.2f} median {np.median(s):.2f} "
              f"max {s.max():.2f}", flush=True)
        if nr:  # per-cell norm restoration: fused change rescaled toward the typical single-source size, ^nr
            srcs = rows if cd4 > 0 else rows[:-1]
            norms = np.array([np.linalg.norm(np.where(np.isfinite(x) & ~own, x, 0), axis=1) for x in srcs])
            have = np.array([np.isfinite(x).any(1) for x in srcs])
            single = (norms * have).sum(0) / np.maximum(have.sum(0), 1)
            r = np.clip(single / (np.linalg.norm(np.where(own, 0, ec), axis=1) + 1e-9), 1, nrmax) ** nr
            print(f"  norm restoration: ratio median {np.median(r):.2f}, max {r.max():.2f}", flush=True)
            ec = np.where(own, ec, ec * r[:, None])
        if gconf:  # per gene: sign agreement across sources |sum_s e_s| / sum_s |e_s| in [0,1], ^gconf, target norm kept
            srcs = rows if cd4 > 0 else rows[:-1]
            num = np.abs(sum(np.nan_to_num(x) for x in srcs))
            den = sum(np.abs(np.nan_to_num(x)) for x in srcs)
            nsrc = sum(np.isfinite(x) for x in srcs)
            agr = np.where(nsrc >= 2, num / np.maximum(den, 1e-12), 0.5)
            w = agr ** gconf
            E0 = np.where(own, 0, ec)
            E1 = E0 * w
            k = np.linalg.norm(E0, axis=1) / (np.linalg.norm(E1, axis=1) + 1e-12)
            print(f"  gene confidence: median agreement {np.median(agr[nsrc >= 2]):.2f}, norm rescale median {np.median(k):.2f}", flush=True)
            ec = np.where(own, ec, E1 * k[:, None])
        if ew:  # shrink per-cell changes of genes weakly expressed in this destination: x sqrt(cpm / (cpm + ew))
            cpm = 1e6 * np.asarray(m) / np.asarray(m).sum()
            f = np.sqrt(cpm / (cpm + ew))
            ec = np.where(own, ec, ec * f[None, :])
        if eb_boost:  # boost per-cell changes of well-expressed genes: x (1 + eb_boost * sqrt(cpm / (cpm + ek)))
            cpm = 1e6 * np.asarray(m) / np.asarray(m).sum()
            f = 1 + eb_boost * np.sqrt(cpm / (cpm + ek))
            ec = np.where(own, ec, ec * f[None, :])
        if pc:  # panel centering of the pooled change: subtract pc x the mean over this panel's targets (own genes excluded)
            E = np.where(own, np.nan, eb)
            mu = np.nanmean(E, axis=0)
            eb = np.where(own, eb, eb - pc * np.nan_to_num(mu)[None, :])
            print(f"  panel centering pc={pc}: mean-change norm {np.linalg.norm(np.nan_to_num(mu)):.3f} vs median target "
                  f"{np.median(np.linalg.norm(np.where(own, 0, eb), axis=1)):.3f}", flush=True)
        if cq:  # DE-consensus on the per-cell change: top-cq |fused| genes whose sources all agree in sign x cboost, rest x csupp
            srcs = rows if cd4 > 0 else rows[:-1]
            pos = sum((np.nan_to_num(x) > 0).astype(int) for x in srcs)
            neg = sum((np.nan_to_num(x) < 0).astype(int) for x in srcs)
            nsrc = sum(np.isfinite(x).astype(int) for x in srcs)
            agree = (nsrc >= cmin) & ((pos == nsrc) | (neg == nsrc))
            E = np.abs(np.where(own, 0, ec))
            thr = np.quantile(E, 1 - cq, axis=1, keepdims=True)
            top = agree & (E >= thr)
            print(f"  DE consensus: boosted genes per target median {np.median(top.sum(1)):.0f}", flush=True)
            ec = np.where(own, ec, ec * np.where(top, cboost, csupp))
        if pcl or pcg:  # context PCA smoothing of the pooled change (and of the per-cell change if pcc)
            eb = pca_smooth(eb, own, pcl, pcn, pcg, pclo, pcspec, pchi)
            if pcc:
                ec = pca_smooth(ec, own, pcc, pcn, pcg, pclo, pcspec, pchi)
        keep_own = own  # the knocked-down gene keeps its full change
        eb = np.where(keep_own, eb, eb * s[:, None])
        if cell:
            ec = np.where(keep_own, ec, ec * s[:, None])
        dc = A.desired_mean(m, ec, space="log2fc", amplitude=ac, clip=3)
        db = A.desired_mean(q, eb, space="bulk_delta", amplitude=ab, clip=3)
        pairs = A.promoter_pairs(targets, genes)
        dc, _ = A.apply_promoter_prior(dc, m, np.asarray(targets), genes, pairs, .15)
        db, _ = A.apply_promoter_prior(db, q, np.asarray(targets), genes, pairs, .15)
        return dc, db
    return profiles


if __name__ == "__main__":
    cmd, rest = sys.argv[1], sys.argv[2:]
    if cmd == "local":
        line, name, rest = rest[0], rest[1], rest[2:]
    else:
        line, name, rest = "2026", rest[0], rest[1:]
    o = {k: float(v) for k, v in (a.split("=") for a in rest)}
    ac, ab = o.pop("ac", 1.0), o.pop("ab", 0.5)
    pool = int(o.pop("pool", 4))
    for src in ("kolf", "hct116", "hek293t", "k562", "h1"):  # source weight overrides, e.g. kolf=1
        if src in o:
            w = o.pop(src)
            A.WEIGHTS = {**A.WEIGHTS, src: w} if w > 0 else {k: v for k, v in A.WEIGHTS.items() if k != src}
    pcs = int(o.pop("pcs", 0))
    if pcs:
        capture_controls(pcs)
    tpow = o.pop("tpow", 1.0)
    if tpow != 1.0:  # Codex's template-variance lever (src/atlas_template_variance.py): templates^power around their mean
        original = A.control_stats

        def stats(*a, **kw):
            m, q, tpl, depths = original(*a, **kw)
            mean = tpl.mean(0)
            ratio = np.divide(tpl, mean, out=np.zeros_like(tpl), where=mean > 0)
            varied = mean * np.power(ratio, tpow)
            return m, q, varied / varied.sum(1, keepdims=True), depths
        A.control_stats = stats
    if os.environ.get("ATLAS_OUT"):  # e.g. Codex's VIP blend cache (data/viperturb/atlas_mix_*); built on atlas_shift_x
        from pathlib import Path
        A.OUT = Path(os.environ["ATLAS_OUT"])
        C.PATH = A.ROOT / "data/atlas_shift_x/cd4_de.npz"
        print(f"source caches from {A.OUT}", flush=True)
    if o.pop("x", 0):
        import context_weights as CW
        CW.use_x()
    A.profiles = make_profiles(**o)
    A.dual_moment_counts = soft_generator()
    print(f"{line}/{name}: agreement allocation {o}, ac={ac}, ab={ab}", flush=True)
    if cmd == "local":
        A.local(line, name, ac=ac, ab=ab, pool=pool)
    else:
        A.build_2026(name, ac=ac, ab=ab, pool=pool)
    fb = A.dual_moment_counts.fallbacks
    print(f"generator fallbacks: {len(fb)} targets, shrink factors {sorted(fb)[:20]}", flush=True)
