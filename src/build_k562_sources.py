"""Alternative K562 sources: vetted controls and batch-normalized effects.

The authors' normalized bulk values are gemgroup Z scores, not counts. We never
submit those as counts. A control-fitted scale maps their deltas back to a count
proxy for the existing transfer pipeline. This mapping is experimental and must
win local comparisons. Original source files are preserved.
"""
import json

import anndata as ad
import numpy as np

from lines import ROOT, CTRL, dedup, weighted_by_target


def main():
    raw = ad.read_h5ad(ROOT / "data/replogle/K562_gwps_raw_bulk_01.h5ad", backed="r")
    norm = ad.read_h5ad(ROOT / "data/replogle/K562_gwps_normalized_bulk_01.h5ad", backed="r")
    assert raw.obs_names.equals(norm.obs_names) and raw.var_names.equals(norm.var_names)
    labels = np.array([name.split("_")[1] for name in raw.obs_names])
    weights = raw.obs.num_cells_filtered.fillna(raw.obs.num_cells_unfiltered).to_numpy(float)
    core = np.flatnonzero((labels == CTRL) & raw.obs.core_control.to_numpy(bool))
    if len(core) < 20:
        raise ValueError("Too few vetted control guides")
    ctrl_counts = np.asarray(raw.X[core], dtype=np.float64)
    ctrl_weights = weights[core]
    ctrl = np.average(ctrl_counts, axis=0, weights=ctrl_weights)
    # Normalize guide depths before estimating the count-space gene scales.
    Y = ctrl_counts * (ctrl.sum() / ctrl_counts.sum(1))[:, None]
    ctrl_normalized = np.average(Y, axis=0, weights=ctrl_weights)
    raw_variance = (ctrl_weights[:, None] * (Y - ctrl_normalized) ** 2).sum(0) / (len(core) - 1)
    Zctrl = np.asarray(norm.X[core], dtype=np.float64)
    valid = np.isfinite(Zctrl)
    finite_weights = valid * ctrl_weights[:, None]
    safe = np.where(valid, Zctrl, 0)
    zmean = (safe * ctrl_weights[:, None]).sum(0) / np.maximum(finite_weights.sum(0), 1)
    zvar = (finite_weights * (safe - zmean) ** 2).sum(0) / np.maximum(valid.sum(0) - 1, 1)
    scale = np.sqrt(raw_variance / np.maximum(zvar, 1e-12))
    scale[(valid.sum(0) < 20) | (zvar < 1e-8)] = 0
    genes = np.array(raw.var.gene_name, dtype=str)

    original = dict(np.load(ROOT / "data/lines/k562.npz"))
    core_genes, core_ctrl = dedup(genes, ctrl)
    assert np.array_equal(core_genes, original["genes"])
    original["ctrl"] = core_ctrl
    original["ctrl_n"] = ctrl_weights.sum()
    np.savez_compressed(ROOT / "data/lines/k562_core.npz", **original)
    del original

    Z = np.asarray(norm.X[:], dtype=np.float32)
    nonfinite = int((~np.isfinite(Z)).sum())
    # Invalid normalizations provide no estimated effect; do not turn infinities
    # into enormous artificial gene changes.
    Z = np.where(np.isfinite(Z), Z, zmean)
    targets, zcounts, n = weighted_by_target(labels, Z, weights)
    counts = np.maximum(ctrl_normalized + (zcounts - zmean) * scale, 0)
    keep = targets != CTRL
    genes, counts, batch_ctrl = dedup(genes, counts[keep], ctrl_normalized)
    np.savez_compressed(ROOT / "data/lines/k562_batch.npz", genes=genes,
                        targets=targets[keep], counts=counts, n=n[keep],
                        ctrl=batch_ctrl, ctrl_n=ctrl_weights.sum())
    summary = dict(source_article="https://api.figshare.com/v2/articles/20029387",
                   raw_file_id=35774443, normalized_file_id=35773217,
                   core_guides=len(core), rejected_control_guides=int((labels == CTRL).sum() - len(core)),
                   genes=len(genes), targets=int(keep.sum()), nonfinite_normalized_entries=nonfinite,
                   zero_scale_genes=int((scale == 0).sum()),
                   mapping="vetted depth-normalized control variance / normalized control variance",
                   status="experimental count proxy; not original observed counts")
    (ROOT / "data/calibration/k562_sources.json").write_text(json.dumps(summary, indent=2))
    raw.file.close()
    norm.file.close()
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
