"""Stream local prediction blocks without retaining an entire cell matrix."""
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def write_prediction_blocks(path, genes, labels, blocks):
    obs = pd.DataFrame({"target": labels}, index=[f"pred_{i}" for i in range(len(labels))])
    ad.AnnData(obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(path)
    rows = 0
    with h5py.File(path, "a") as f:
        matrix = None
        for block in blocks:
            block = sp.csr_matrix(block, dtype=np.float32)
            if block.shape[1] != len(genes):
                raise ValueError("Prediction block has wrong gene dimension")
            rows += block.shape[0]
            if rows > len(labels):
                raise ValueError("More prediction rows than labels")
            if matrix is None:
                ad.io.write_elem(f, "X", block)
                matrix = ad.io.sparse_dataset(f["X"])
            else:
                matrix.append(block)
        if rows != len(labels):
            raise ValueError(f"Prediction has {rows} rows for {len(labels)} labels")
