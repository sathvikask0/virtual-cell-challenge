"""State pilot on Modal: retrain Arc State (ST, HVG space) with ESM2 perturbation features, H1 held out zero-shot.

Steps (from src/, each is a separate `modal run` so costs stay visible):
  modal volume put vcc-state ../data/state_pilot /raw              # local subsampled jurkat/hepg2/h1 + esm2.pt
  modal run modal_state.py::fetch --name k562gw                    # CPU: download 66 GB, keep <=40 cells/target
  modal run modal_state.py::fetch --name rpe1
  modal run modal_state.py::build                                  # CPU: normalise, 2,000 shared HVGs, TOML
  modal run modal_state.py::train --steps 20000                    # GPU L40S
  modal run modal_state.py::predict                                # GPU: H1 predictions -> /runs/pilot1/...
  modal volume get vcc-state /runs/pilot1/eval_last.ckpt/adata_pred.h5ad ../data/state_pilot/h1_pred.h5ad
"""
import modal

app = modal.App("vcc-state")
vol = modal.Volume.from_name("vcc-state", create_if_missing=True)
cpu_img = modal.Image.debian_slim(python_version="3.12").pip_install(
    "anndata==0.11.4", "h5py==3.16.0", "numpy<2.3", "pandas<3", "scipy", "scanpy==1.10.4", "typing_extensions>=4.12", "toml")
gpu_img = modal.Image.debian_slim(python_version="3.12").pip_install("arc-state==0.11.1", "h5py==3.16.0")

URLS = {"k562gw": "https://ndownloader.figshare.com/files/35775507",   # K562_gwps_raw_singlecell_01.h5ad
        "rpe1": "https://ndownloader.figshare.com/files/35775606"}     # rpe1_raw_singlecell_01.h5ad
PER_TARGET, N_CTRL, BLOCK = 40, 20000, 50000


@app.function(image=cpu_img, volumes={"/vol": vol}, cpu=4, memory=32768, timeout=4 * 3600,
              ephemeral_disk=512 * 1024)
def fetch(name: str):
    """Download a Replogle raw single-cell h5ad to local disk, stream it in row blocks, keep a subsample."""
    import shutil
    import urllib.request

    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from anndata.io import read_elem

    path = f"/tmp/{name}.h5ad"
    for attempt in range(5):
        try:
            with urllib.request.urlopen(URLS[name], timeout=600) as r, open(path, "wb") as f:
                shutil.copyfileobj(r, f, length=64 << 20)
            break
        except Exception as e:
            print(f"  download attempt {attempt + 1} failed: {e}", flush=True)
    else:
        raise RuntimeError(f"{name}: download failed")
    with h5py.File(path, "r") as f:
        obs, var = read_elem(f["obs"]), read_elem(f["var"])
        lab = obs["gene"].astype(str).to_numpy()
        rng = np.random.default_rng(0)
        keep = []
        for t, idx in pd.Series(np.arange(len(lab))).groupby(lab):
            n = N_CTRL if t == "non-targeting" else PER_TARGET
            idx = idx.to_numpy()
            keep.append(idx if len(idx) <= n else rng.choice(idx, n, replace=False))
        keep = np.sort(np.concatenate(keep))
        X = f["X"]
        rows = []
        for s in range(0, len(lab), BLOCK):
            sel = keep[(keep >= s) & (keep < s + BLOCK)]
            if len(sel) == 0:
                continue
            if not isinstance(X, h5py.Dataset):
                raise ValueError(f"{name}: expected a dense X (Replogle raw single-cell files are dense)")
            rows.append(sp.csr_matrix(X[s:min(s + BLOCK, len(lab))][sel - s].astype(np.float32)))
            print(f"  {name}: rows {s:,}/{len(lab):,}", flush=True)
    M = sp.vstack(rows).tocsr()
    genes = var["gene_name"].astype(str).to_numpy() if "gene_name" in var else var.index.astype(str).to_numpy()
    o = pd.DataFrame({"gene": np.asarray(lab[keep], dtype=object), "cell_type": name,
                      "gem_group": np.asarray(obs["gem_group"].astype(str).to_numpy()[keep], dtype=object)},
                     index=pd.Index([str(i) for i in keep], dtype=object))
    a = ad.AnnData(M, obs=o, var=pd.DataFrame(index=pd.Index(np.asarray(genes, dtype=object), dtype=object)))
    a.var_names_make_unique()
    a.write_h5ad(f"/vol/raw/{name}.h5ad", compression="gzip")
    vol.commit()
    print(f"{name}: {a.n_obs:,} cells, {a.n_vars:,} genes, {o.gene.nunique():,} perturbations")


@app.function(image=cpu_img, volumes={"/vol": vol}, cpu=4, memory=65536, timeout=3 * 3600)
def build(n_hvg: int = 2000, test: str = "h1"):
    """Normalise (CP10k log1p), pick shared HVGs, write State inputs (X_hvg) and the TOML (test line zero-shot)."""
    import anndata as ad
    import numpy as np
    import scanpy as sc
    import toml

    names = ["k562gw", "rpe1", "jurkat", "hepg2", "h1"]
    import json
    import os
    A = {n: ad.read_h5ad(f"/vol/raw/{n}.h5ad") for n in names}
    for n in ("k562gw", "rpe1"):  # legacy categorical var in the Replogle files read back as codes: restore symbols
        gfile = f"/vol/raw/{n}_genes.json"
        if os.path.exists(gfile):
            g = json.load(open(gfile))
            assert len(g) == A[n].n_vars, (n, len(g), A[n].n_vars)
            A[n].var_names = g
            A[n].var_names_make_unique()
            print(f"  {n}: gene names restored, e.g. {g[:3]}", flush=True)
    for n, a in A.items():
        print(f"  {n}: {a.shape}, obs columns {list(a.obs.columns)}", flush=True)
        if "gene" not in a.obs.columns:
            raise ValueError(f"{n}: no obs 'gene' column")
    common = sorted(set.intersection(*(set(a.var_names) for a in A.values())))
    print(f"shared genes {len(common)}", flush=True)
    pooled = []
    for n, a in A.items():
        a = a[:, common].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        A[n] = a
        c = a[a.obs.gene == "non-targeting"]
        pooled.append(c[np.random.default_rng(0).choice(c.n_obs, min(5000, c.n_obs), replace=False)])
    P = ad.concat(pooled, label="dataset", keys=names)
    sc.pp.highly_variable_genes(P, n_top_genes=n_hvg, batch_key="dataset")
    hvg = P.var_names[P.var.highly_variable].tolist()
    import os
    os.makedirs("/vol/proc", exist_ok=True)
    for n, a in A.items():
        b = a[:, hvg].copy()
        b.obsm["X_hvg"] = np.asarray(b.X.toarray() if hasattr(b.X, "toarray") else b.X, np.float32)
        # State's loader reads obs/_index and var/_index as plain string datasets: no nullable strings
        import pandas as pd
        b.obs_names = pd.Index(np.asarray(b.obs_names, dtype=object), dtype=object)
        b.var_names = pd.Index(np.asarray(b.var_names, dtype=object), dtype=object)
        for c in b.obs.columns:
            b.obs[c] = pd.Categorical(np.asarray(b.obs[c], dtype=object))
        b.write_h5ad(f"/vol/proc/{n}.h5ad")
        print(f"  {n}: {b.n_obs:,} cells", flush=True)
    cfg = {"datasets": {n: f"/vol/proc/{n}.h5ad" for n in names},
           "training": {n: "train" for n in names},
           "zeroshot": {f"{test}.{test}": "test", "hepg2.hepg2": "val"}}
    with open("/vol/proc/pilot.toml", "w") as f:
        toml.dump(cfg, f)
    with open("/vol/proc/hvg.txt", "w") as f:
        f.write("\n".join(hvg))
    vol.commit()


@app.function(image=gpu_img, volumes={"/vol": vol}, gpu="L40S", cpu=8, memory=65536, timeout=8 * 3600)
def train(steps: int = 20000, name: str = "pilot1"):
    import subprocess
    subprocess.run([
        "state", "tx", "train",
        "data.kwargs.toml_config_path=/vol/proc/pilot.toml", "data.kwargs.embed_key=X_hvg",
        "data.kwargs.output_space=gene", "data.kwargs.pert_col=gene", "data.kwargs.cell_type_key=cell_type",
        "data.kwargs.batch_col=gem_group", "data.kwargs.control_pert=non-targeting",
        "data.kwargs.perturbation_features_file=/vol/raw/esm2.pt", "data.kwargs.num_workers=6",
        f"training.max_steps={steps}", "training.val_freq=2000", "training.ckpt_every_n_steps=5000",
        "model=state_sm", "model.kwargs.hidden_dim=328", "model.kwargs.cell_set_len=64",
        "use_wandb=false", "output_dir=/vol/runs", f"name={name}"], check=True)
    vol.commit()


@app.function(image=gpu_img, volumes={"/vol": vol}, gpu="L40S", cpu=8, memory=65536, timeout=2 * 3600)
def predict(name: str = "pilot1", ckpt: str = "last.ckpt"):
    import subprocess
    subprocess.run(["state", "tx", "predict", "--output-dir", f"/vol/runs/{name}", "--checkpoint", ckpt,
                    "--profile", "anndata"], check=True)
    vol.commit()
