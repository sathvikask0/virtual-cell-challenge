"""Control-fitted low-rank Poisson-lognormal sampler with fixed marginal moments.

The independent and factor variants share rate mean=1 and variance=phi.
Factors are learned from input controls only; no reference perturbations used.
"""
import numpy as np
import scipy.sparse as sp
from sklearn.utils.extmath import randomized_svd


def fit_factors(X, phi, rank=32, max_cells=2000, seed=0):
    X = sp.csr_matrix(X, dtype=np.float64)
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(X.shape[0], min(max_cells, X.shape[0]), replace=False))
    X = X[rows]
    totals = np.asarray(X.sum(1)).ravel()
    Y = sp.diags(np.median(totals) / np.maximum(totals, 1)) @ X
    mean = np.asarray(Y.mean(0)).ravel()
    var = np.asarray(Y.multiply(Y).mean(0)).ravel() - mean ** 2
    keep = np.flatnonzero((mean > .2) & (var > 0) & (phi > 0))
    rank = min(rank, len(keep), len(rows) - 1)
    if rank < 1:
        return np.zeros((X.shape[1], 0))
    Z = Y[:, keep].toarray().astype(np.float32)
    Z -= mean[keep].astype(np.float32)
    Z /= np.sqrt(var[keep]).astype(np.float32)
    _, singular, components = randomized_svd(Z, n_components=rank, n_iter=3, random_state=seed)
    loading = components.T * (singular / np.sqrt(len(rows) - 1))
    loading *= (np.sqrt(var[keep]) / mean[keep])[:, None]
    # Bound every gene's shared variance by its fitted biological rate variance.
    sigma2 = np.log1p(phi[keep])
    shared = np.square(loading).sum(1)
    loading *= np.sqrt(np.minimum(1., sigma2 / np.maximum(shared, 1e-12)))[:, None]
    result = np.zeros((X.shape[1], rank), dtype=np.float32)
    result[keep] = loading
    return result


def sample_lognormal(share, totals, phi, ltr, rng, n=400, factors=None, strength=1.):
    sigma2 = np.log1p(phi)
    latent = np.zeros((n, len(share)))
    shared = np.zeros(len(share))
    if factors is not None and factors.shape[1]:
        loading = factors * strength
        shared = np.square(loading).sum(1)
        if np.any(shared > sigma2 + 1e-6):
            raise ValueError("Factor variance exceeds marginal variance")
        latent = rng.normal(size=(n, loading.shape[1])) @ loading.T
    latent += rng.normal(size=latent.shape) * np.sqrt(np.maximum(sigma2 - shared, 0))
    rate = np.exp(latent - sigma2 / 2)
    lib = rng.choice(totals, n) * np.exp(ltr)
    return sp.csr_matrix(rng.poisson(lib[:, None] * share * rate).astype(np.float32))
