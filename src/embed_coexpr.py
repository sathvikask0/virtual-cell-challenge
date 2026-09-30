"""Gene embeddings from co-expression in control cells.

Genes that rise and fall together across the 38k control cells get nearby vectors.
Only control (non-targeting) cells are used, so no knockdown answers leak in.

Steps:
  1. each control cell: scale to 10k total counts, log1p
  2. z-score each gene across cells
  3. truncated SVD of the cells x genes matrix; gene vector = right singular vector * singular value
     (so dot products between gene vectors approximate their correlation)

Output: data/embeddings/coexpr.npz with genes (18080,), emb (18080, DIM), sv (DIM,)
"""
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.utils.extmath import randomized_svd

ROOT = Path(__file__).resolve().parent.parent
CHUNK = 20_000
DIM = 128


def main():
    a = ad.read_h5ad(ROOT / "data/vcc/train/adata_Training.h5ad", backed="r")
    is_ctrl = (a.obs["target_gene"] == "non-targeting").to_numpy()

    parts = []
    for start in range(0, a.n_obs, CHUNK):
        stop = min(start + CHUNK, a.n_obs)
        keep = is_ctrl[start:stop]
        if not keep.any():
            continue
        X = sp.csr_matrix(a.X[start:stop])[keep].astype(np.float32)
        X = sp.diags(1e4 / np.maximum(X.sum(1).A1, 1)) @ X
        X.data = np.log1p(X.data)
        parts.append(X.toarray())
        print(f"  read {stop:,}/{a.n_obs:,} cells", flush=True)
    Z = np.concatenate(parts)
    del parts
    print(f"  control matrix: {Z.shape[0]:,} cells x {Z.shape[1]:,} genes")

    Z -= Z.mean(0)
    sd = Z.std(0)
    Z /= np.where(sd > 0, sd, 1)

    _, S, Vt = randomized_svd(Z, n_components=DIM, random_state=0)
    emb = (Vt.T * S / np.sqrt(Z.shape[0])).astype(np.float32)

    out = ROOT / "data/embeddings/coexpr.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, genes=np.array(a.var_names, dtype=str), emb=emb, sv=S)
    print(f"  saved {out} ({emb.shape[0]:,} genes x {DIM} dims)")


if __name__ == "__main__":
    main()
