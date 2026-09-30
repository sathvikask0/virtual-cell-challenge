"""Sum raw counts per knocked-down gene (and per gene x batch).

For each group we keep three things, from which everything else follows:
  sum_counts  (18080,)  raw counts of each measured gene, summed over the group's cells
  sum_total   scalar    total counts, summed over the group's cells
  n           scalar    number of cells

  share       = sum_counts / sum_total   what fraction of the group's RNA each gene is
  mean_total  = sum_total / n            average total counts per cell
  mean_counts = sum_counts / n           = share * mean_total, the counts-space answer

Output: data/pseudobulk/{split}.npz with genes, targets, sum_counts, sum_total, n,
and the same per (target, batch) as tb_target, tb_batch, tb_sum_counts, tb_sum_total, tb_n.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
FILES = {
    "train": "train/adata_Training.h5ad",
    "val": "validation/adata_Validation.h5ad",
    "test": "test/adata_Test.h5ad",
}
CHUNK = 20_000


def group_sums(X, codes, n_groups):
    """Sum rows of X by integer group code (via a sparse indicator matmul)."""
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                      shape=(n_groups, len(codes)))
    return np.asarray((G @ X).todense())


def run(split):
    a = ad.read_h5ad(ROOT / "data/vcc" / FILES[split], backed="r")
    obs = a.obs
    t_cat = obs["target_gene"].astype("category")
    tb_cat = (obs["target_gene"].astype(str) + "|" + obs["batch"].astype(str)).astype("category")
    t_codes, tb_codes = t_cat.cat.codes.to_numpy(), tb_cat.cat.codes.to_numpy()
    nT, nTB, nG = len(t_cat.cat.categories), len(tb_cat.cat.categories), a.n_vars

    t_sum, tb_sum = np.zeros((nT, nG)), np.zeros((nTB, nG))
    for start in range(0, a.n_obs, CHUNK):
        stop = min(start + CHUNK, a.n_obs)
        X = sp.csr_matrix(a.X[start:stop], dtype=np.float64)
        t_sum += group_sums(X, t_codes[start:stop], nT)
        tb_sum += group_sums(X, tb_codes[start:stop], nTB)
        print(f"  {split}: {stop:,}/{a.n_obs:,} cells", flush=True)

    tb_labels = np.array(tb_cat.cat.categories, dtype=str)
    out = ROOT / "data/pseudobulk" / f"{split}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        genes=np.array(a.var_names, dtype=str),
        targets=np.array(t_cat.cat.categories, dtype=str),
        sum_counts=t_sum, sum_total=t_sum.sum(1), n=np.bincount(t_codes, minlength=nT),
        tb_target=np.array([s.split("|")[0] for s in tb_labels]),
        tb_batch=np.array([s.split("|")[1] for s in tb_labels]),
        tb_sum_counts=tb_sum, tb_sum_total=tb_sum.sum(1),
        tb_n=np.bincount(tb_codes, minlength=nTB),
    )
    print(f"  saved {out} ({nT} targets, {nTB} target x batch groups)")


if __name__ == "__main__":
    for split in sys.argv[1:] or FILES:
        run(split)
