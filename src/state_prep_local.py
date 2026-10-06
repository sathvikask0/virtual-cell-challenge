"""State pilot, local half: subsample Nadig Jurkat/HepG2 and H1 into small raw-count h5ads for Modal.

Each output has raw counts in X and obs columns gene (target, 'non-targeting' for controls), cell_type, gem_group.
Up to PER_TARGET cells per target, up to N_CTRL controls. H1 is the held-out (zero-shot) test line: its 150
local-eval targets plus controls; State never trains on it.
Usage (from src/): python state_prep_local.py   -> data/state_pilot/{jurkat,hepg2,h1}.h5ad
"""
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import atlas_shift as A

OUT = A.ROOT / "data/state_pilot"
PER_TARGET, N_CTRL = 40, 20000
rng = np.random.default_rng(0)


def pick(labels, ctrl="non-targeting"):
    keep = []
    for t, idx in pd.Series(np.arange(len(labels))).groupby(labels):
        n = N_CTRL if t == ctrl else PER_TARGET
        keep.append(idx.to_numpy() if len(idx) <= n else rng.choice(idx.to_numpy(), n, replace=False))
    return np.sort(np.concatenate(keep))


def write(name, X, genes, gene, cell_type, batch):
    obs = pd.DataFrame({"gene": np.asarray(gene, str), "cell_type": cell_type,
                        "gem_group": np.asarray(batch, str)})
    a = ad.AnnData(sp.csr_matrix(X, dtype=np.float32), obs=obs, var=pd.DataFrame(index=np.asarray(genes, str)))
    a.var_names_make_unique()
    a.write_h5ad(OUT / f"{name}.h5ad", compression="gzip")
    print(f"{name}: {a.n_obs:,} cells, {a.n_vars:,} genes, {a.obs.gene.nunique():,} perturbations", flush=True)


def nadig(line):
    a = ad.read_h5ad(A.ROOT / f"data/nadig/{line}_raw_singlecell.h5ad", backed="r")
    lab = a.obs["gene"].astype(str).to_numpy()
    idx = pick(lab)
    X = a.X[idx]
    write(line, X, a.var["gene_name"] if "gene_name" in a.var else a.var_names, lab[idx], line,
          a.obs["gem_group"].astype(str).to_numpy()[idx])


def h1():
    r = ad.read_h5ad(A.ROOT / "data/local_eval/h1/real.h5ad")
    c = ad.read_h5ad(A.ROOT / "data/local_eval/h1/ctrl_input.h5ad")
    lab = r.obs["target"].astype(str).to_numpy()
    pert = np.flatnonzero(lab != "non-targeting")
    pert = pick(lab[pert]) if False else pert[pick(lab[pert])]
    cidx = rng.choice(c.n_obs, min(N_CTRL, c.n_obs), replace=False)
    assert list(c.var_names) == list(r.var_names)
    X = sp.vstack([sp.csr_matrix(r.X[pert]), sp.csr_matrix(c.X[cidx])])
    gene = np.concatenate([lab[pert], np.full(len(cidx), "non-targeting")])
    write("h1", X, r.var_names, gene, "h1", np.zeros(len(gene), int))


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    h1()
    for line in ("jurkat", "hepg2"):
        nadig(line)
