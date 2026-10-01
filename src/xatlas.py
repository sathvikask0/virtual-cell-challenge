"""X-Atlas/Orion (Xaira 2025): genome-wide CRISPRi in HCT116 and HEK293T, into the lines format.

The release is per-cell only (~127 GB of parquet on Hugging Face, every batch mixes all targets),
so this streams it: download one batch (~0.4 GB) while summing the previous one, delete it, repeat.
Only the 2026 genes are kept. Progress is saved every SAVE_EVERY batches, so a rerun resumes.

Output (same format as src/lines.py and src/noise.py):
  data/lines/{hct116,hek293t}.npz        genes, targets, counts (mean per cell), n, ctrl, ctrl_n
  data/lines/{hct116,hek293t}_noise.npz  genes, var (per control cell)

Usage: python src/xatlas.py [HCT116|HEK293T ...]   (default both)
"""
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/xatlas"
REPO = "https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion"
CTRL = "Non-Targeting"
SAVE_EVERY = 10


def batches(line):
    with urllib.request.urlopen(f"https://huggingface.co/api/datasets/Xaira-Therapeutics/"
                                f"X-Atlas-Orion/tree/main/data") as r:
        files = [f["path"] for f in json.load(r) if f["path"].startswith(f"data/{line}_Batch")]
    return sorted(files, key=lambda p: int(p.split("Batch")[1].split(".")[0]))


def fetch(path, dest):
    return subprocess.Popen(["curl", "-sSL", "--retry", "5", "--retry-delay", "10",
                             "-o", str(dest), f"{REPO}/resolve/main/{path}"])


def add_batch(path, col_of_token, G, acc):
    """Add one parquet batch to the running sums in acc."""
    pf = pq.ParquetFile(path)
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=["gene_token_id", "gene_expression", "gene_target",
                                           "pass_guide_filter"])
        keep = np.asarray(t["pass_guide_filter"]) == 1
        targets = np.asarray(t["gene_target"].to_pylist(), dtype=object)[keep]
        tok, val = t["gene_token_id"].filter(keep), t["gene_expression"].filter(keep)
        offsets = np.asarray(tok.combine_chunks().offsets)
        cols = col_of_token[np.asarray(pc.list_flatten(tok))]
        vals = np.asarray(pc.list_flatten(val), dtype=np.float32)
        rows = np.repeat(np.arange(len(targets)), np.diff(offsets))
        m = cols >= 0
        X = sp.csr_matrix((vals[m], (rows[m], cols[m])), shape=(len(targets), G))
        for tg in np.unique(targets):
            if tg not in acc["row"]:
                acc["row"][tg] = len(acc["row"])
        ridx = np.array([acc["row"][tg] for tg in targets])
        A = sp.csr_matrix((np.ones(len(ridx), np.float32), (ridx, np.arange(len(ridx)))),
                          shape=(len(acc["row"]), len(ridx)))
        S = (A @ X).tocoo()
        if acc["sum"].shape[0] < len(acc["row"]):  # grow the target x gene sum matrix
            grown = np.zeros((len(acc["row"]) + 2000, G), np.float32)
            grown[:acc["sum"].shape[0]] = acc["sum"]
            acc["sum"] = grown
        np.add.at(acc["sum"], (S.row, S.col), S.data)
        if len(acc["n"]) < len(acc["row"]):
            acc["n"] = np.concatenate([acc["n"], np.zeros(len(acc["row"]) - len(acc["n"]))])
        np.add.at(acc["n"], ridx, 1)
        c = targets == CTRL
        if c.any():
            Xc = X[c]
            acc["ctrl_sq"] += np.asarray(Xc.multiply(Xc).sum(0), dtype=np.float64).ravel()


def run(line):
    genes = pd.read_csv(ROOT / "data/vcc/controls/gene_names.csv")["gene_name"].to_numpy(str)
    meta = pd.read_parquet(OUT / "gene_metadata.parquet")
    gidx = {g: i for i, g in enumerate(genes)}
    col_of_token = np.full(meta["gene_token_id"].max() + 1, -1)
    for tok, name in zip(meta["gene_token_id"], meta["gene_name"]):
        if name in gidx and col_of_token[tok] < 0:
            col_of_token[tok] = gidx[name]
    G = len(genes)
    part = OUT / f"partial_{line}.npz"
    if part.exists():
        p = np.load(part, allow_pickle=True)
        acc = {"sum": p["sum"], "n": p["n"], "ctrl_sq": p["ctrl_sq"],
               "row": {t: i for i, t in enumerate(p["targets"])}}
        done = set(p["done"])
    else:
        acc = {"sum": np.zeros((0, G), np.float32), "n": np.zeros(0), "ctrl_sq": np.zeros(G),
               "row": {}}
        done = set()
    todo = [b for b in batches(line) if b not in done]
    print(f"{line}: {len(done)} batches done, {len(todo)} to go", flush=True)

    def save():
        targets = np.array(sorted(acc["row"], key=acc["row"].get), dtype=object)
        np.savez(part, sum=acc["sum"], n=acc["n"], ctrl_sq=acc["ctrl_sq"], targets=targets,
                 done=np.array(sorted(done), dtype=object))

    nxt = None
    for i, b in enumerate(todo):
        dest = OUT / f"tmp_{line}_{i % 2}.parquet"
        proc = nxt or fetch(b, dest)
        if proc.wait() != 0:
            raise RuntimeError(f"download failed: {b}")
        if i + 1 < len(todo):
            nxt = fetch(todo[i + 1], OUT / f"tmp_{line}_{(i + 1) % 2}.parquet")
        add_batch(dest, col_of_token, G, acc)
        dest.unlink()
        done.add(b)
        if (i + 1) % SAVE_EVERY == 0:
            save()
        print(f"  {line} {len(done)}: {b.split('/')[-1]}, {len(acc['row']):,} targets, "
              f"{int(acc['n'].sum()):,} cells", flush=True)
    save()

    T = len(acc["row"])
    targets = np.array(sorted(acc["row"], key=acc["row"].get), dtype=str)
    S, n = acc["sum"][:T].astype(np.float64), acc["n"][:T]
    ci = acc["row"][CTRL]
    ctrl, ctrl_n = S[ci] / n[ci], n[ci]
    keep = np.arange(T) != ci
    name = line.lower()
    np.savez(ROOT / "data/lines" / f"{name}.npz", genes=genes, targets=targets[keep],
             counts=S[keep] / n[keep, None], n=n[keep], ctrl=ctrl, ctrl_n=ctrl_n)
    var = acc["ctrl_sq"] / ctrl_n - ctrl ** 2
    np.savez(ROOT / "data/lines" / f"{name}_noise.npz", genes=genes, var=var)
    print(f"{name}: {keep.sum():,} targets, median {np.median(n[keep]):.0f} cells, "
          f"{int(ctrl_n):,} controls", flush=True)


if __name__ == "__main__":
    for line in sys.argv[1:] or ["HCT116", "HEK293T"]:
        run(line)
