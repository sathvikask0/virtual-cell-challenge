"""Depth-aware moment fit for normalized-count gamma/Poisson sampling."""
import numpy as np
import scipy.sparse as sp


def depth_aware_stats(X):
    X = sp.csr_matrix(X, dtype=np.float64)
    totals = np.asarray(X.sum(1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("Controls must have positive total counts")
    share = np.asarray(X.sum(0)).ravel() / totals.sum()
    depth = np.median(totals)
    Y = sp.diags(depth / totals) @ X
    m = np.asarray(Y.mean(0)).ravel()
    v = np.asarray(Y.multiply(Y).mean(0)).ravel() - m ** 2
    # Conditional Poisson variance after scaling varies with each cell's depth:
    # E[var(Y_g | L)] = mean(Y_g) * median(L) * E[1/L].
    shot_noise = m * depth * np.mean(1 / totals)
    phi = np.clip((v - shot_noise) / np.maximum(m, 1e-12) ** 2, 0, 10)
    phi[m == 0] = 0
    return share, totals, phi
