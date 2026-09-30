"""Independent transfer calibration; leaves the original submission unchanged."""
import numpy as np

from predict_2026 import EPS, target_change


def calibrated_change(target, gidx, src, share, alpha=.5, total_alpha=.5,
                      own_alpha=1., mode="transfer", expression_gate=0.):
    # Recover the original unshrunk profile, then independently scale its parts.
    lfc, total = target_change(target, gidx, src, 1.)
    own = gidx.get(target)
    own_drop = lfc[own] if own is not None else 0.
    if mode == "control":
        return np.zeros_like(lfc), 0.
    if mode == "target-only":
        lfc = np.zeros_like(lfc)
        total = 0.
    elif mode != "transfer":
        raise ValueError(mode)
    gate = 1.
    if expression_gate > 0 and own is not None:
        # Threshold is normalized expression per 10,000 counts, not raw depth.
        expression = share[own] * 10000
        gate = expression / (expression + expression_gate)
    lfc = alpha * gate * lfc
    if own is not None:
        lfc[own] = own_alpha * own_drop
    return lfc, total_alpha * gate * total


def mean_prediction(target, gidx, src, share, total, **options):
    lfc, ltr = calibrated_change(target, gidx, src, share, **options)
    s = np.maximum((share + EPS) * np.exp(lfc) - EPS, 0)
    return s / s.sum() * total * np.exp(ltr)
