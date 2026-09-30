"""Small learned correction to copied logbulk effects; exploratory outer-line CV.

PCA and ridge are fit exclusively on OTHER lines, including their cross-line
pairs for identical targets. Features combine the source effect with differences
in control expression. No hidden-line perturbation labels enter fitting.

This is a mean-expression proxy, not the full official scorer. All measured genes
are used here, without the official per-cell filter or sampling correction.
Run: python src/learned_transfer.py --lines h1 hepg2 jurkat
"""
import argparse
import json

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from calibration_sweep import ROOT
from cross_line import LABELED, fit, load_lines
from calibration import mean_prediction
from cell_eval2.metrics.discrimination import discrimination_score


def logbulk(counts):
    return np.log1p(50000 * counts / counts.sum(axis=-1, keepdims=True))


def features(source_coeff, context_delta, target_expression):
    interaction = (source_coeff[:, :, None] * context_delta[:, None, :]).reshape(len(source_coeff), -1)
    return np.column_stack([source_coeff, context_delta, interaction, target_expression])


def train_correction(train, genes, rank=32, ridge_alpha=100., use_context=True):
    gidx = {g: i for i, g in enumerate(genes)}
    controls = {n: logbulk(d["ctrl"]) for n, d in train.items()}
    effects = {n: logbulk(d["counts"]) - controls[n] for n, d in train.items()}
    rng = np.random.default_rng(0)
    basis_rows = []
    for n, d in train.items():
        idx = rng.choice(len(d["targets"]), min(500, len(d["targets"])), replace=False)
        rows = effects[n][idx].copy()
        for j, i in enumerate(idx):
            own = gidx.get(d["targets"][i])
            if own is not None:
                rows[j, own] = 0
        basis_rows.append(rows)
    basis = PCA(n_components=rank, svd_solver="randomized", random_state=0)
    basis.fit(np.vstack(basis_rows))
    control_basis = PCA(n_components=min(3, len(train) - 1), svd_solver="full")
    control_basis.fit(np.array(list(controls.values())))
    control_coeff = {n: control_basis.transform(c[None])[0] for n, c in controls.items()}
    # Transform uncentered deltas: correction is a residual with a meaningful zero.
    xs, ys = [], []
    for source, src in train.items():
        for dest, dst in train.items():
            if source == dest:
                continue
            targets = sorted(set(src["targets"]) & set(dst["targets"]))
            if not targets:
                continue
            si = np.array([src["tidx"][t] for t in targets])
            di = np.array([dst["tidx"][t] for t in targets])
            context_delta = np.repeat((control_coeff[dest] - control_coeff[source])[None], len(si), axis=0)
            if not use_context:
                context_delta = np.zeros((len(si), 0))
            expression = np.array([[controls[source][gidx[t]], controls[dest][gidx[t]]]
                                   if t in gidx else [0., 0.] for t in targets])
            source_rows = effects[source][si].copy()
            residual = effects[dest][di] - source_rows
            for j, target in enumerate(targets):
                own = gidx.get(target)
                if own is not None:
                    source_rows[j, own] = 0
                    residual[j, own] = 0
            xs.append(features(source_rows @ basis.components_.T, context_delta, expression))
            ys.append(residual @ basis.components_.T)
    if not xs:
        raise ValueError("No shared training targets for learning cross-line corrections")
    X, y = np.vstack(xs), np.vstack(ys)
    scaler = StandardScaler().fit(X)
    ridge = Ridge(alpha=ridge_alpha).fit(scaler.transform(X), y)
    return dict(basis=basis, controls=controls, effects=effects, control_basis=control_basis,
                control_coeff=control_coeff, scaler=scaler, ridge=ridge, gidx=gidx,
                train_pairs=len(X), feature_count=X.shape[1], use_context=use_context)


def predict_correction(model, train, targets, dest_ctrl):
    dest_control = logbulk(dest_ctrl)
    dest_coeff = model["control_basis"].transform(dest_control[None])[0]
    base_rows, correction_rows = [], []
    for target in targets:
        bases, corrections = [], []
        own = model["gidx"].get(target)
        for source, d in train.items():
            if target not in d["tidx"]:
                continue
            base = model["effects"][source][d["tidx"][target]].copy()
            input_effect = base.copy()
            if own is not None:
                input_effect[own] = 0
            coeff = input_effect @ model["basis"].components_.T
            context = dest_coeff - model["control_coeff"][source]
            if not model["use_context"]:
                context = np.zeros(0)
            expression = [model["controls"][source][own], dest_control[own]] if own is not None else [0., 0.]
            X = features(coeff[None], context[None], np.array(expression)[None])
            residual = model["ridge"].predict(model["scaler"].transform(X))[0] @ model["basis"].components_
            if own is not None:
                residual[own] = 0
            bases.append(base)
            corrections.append(residual)
        if not bases:
            raise ValueError(f"Target has no source: {target}")
        base_rows.append(np.mean(bases, axis=0))
        correction_rows.append(np.mean(corrections, axis=0))
    return np.array(base_rows), np.array(correction_rows)


def evaluate(pred, true, ctrl, targets, genes):
    gidx = {g: i for i, g in enumerate(genes)}
    err = np.square(pred - true)
    control_err = np.square(true - ctrl)
    for row, target in enumerate(targets):
        own = gidx.get(target)
        if own is not None:
            err[row, own] = control_err[row, own] = 0
    labels = np.array([*targets, "non-targeting"])
    pred_bulk = (labels, np.vstack([pred, ctrl]))
    real_bulk = (labels, np.vstack([true, ctrl]))
    pds = discrimination_score(pred_bulk=pred_bulk, real_bulk=real_bulk, genes=genes,
                               distance="cosine", rank_denominator="n-1", exclusion_scope="panel",
                               control_source="real", tie_policy="midrank")
    return dict(relative_logbulk_mse=float(err.sum() / control_err.sum()),
                pds=float(np.mean(list(pds.values()))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lines", nargs="+", choices=LABELED, default=["h1", "hepg2", "jurkat"])
    parser.add_argument("--rank", type=int, default=32)
    parser.add_argument("--ridge-alpha", type=float, default=100.)
    parser.add_argument("--no-context", action="store_true",
                        help="Ablate whole-control PCA and its effect interactions")
    parser.add_argument("--name", default="learned_transfer")
    args = parser.parse_args()
    lines, genes = load_lines(LABELED)
    results = []
    for held in args.lines:
        train = {n: d for n, d in lines.items() if n != held}
        test = lines[held]
        assert held not in train
        model = train_correction(train, genes, args.rank, args.ridge_alpha, not args.no_context)
        available = set().union(*(set(d["targets"]) for d in train.values()))
        targets = sorted(set(test["targets"]) & available)
        rng = np.random.default_rng(0)
        targets = sorted(rng.choice(targets, min(150, len(targets)), replace=False))
        true = logbulk(test["counts"][[test["tidx"][t] for t in targets]])
        ctrl = logbulk(test["ctrl"])
        base, correction = predict_correction(model, train, targets, test["ctrl"])
        for strength in [0., .25, .5, 1.]:
            counts = np.expm1(np.maximum(ctrl + base + strength * correction, 0))
            pred = logbulk(counts)
            scores = evaluate(pred, true, ctrl, targets, genes)
            row = dict(line=held, method="learned_residual", strength=strength,
                       n_targets=len(targets), train_pairs=model["train_pairs"],
                       feature_count=model["feature_count"], **scores)
            results.append(row)
            print(row, flush=True)
        # Compare against the actual EPS-based copying scheme, not just our new base.
        m = fit(train)
        src = ({t: np.mean([x for x, _ in v], axis=0) for t, v in m["per_target"].items()},
               {t: np.ones(len(genes)) for t in m["per_target"]},
               {t: [y for _, y in v] for t, v in m["per_target"].items()},
               m["mean_lfc"], m["mean_ltr"], m["own_drop"])
        pred = logbulk(np.array([mean_prediction(t, m["gidx"], src, test["ctrl_share"],
                                                test["ctrl"].sum(), alpha=1., total_alpha=1.)
                                 for t in targets]))
        scores = evaluate(pred, true, ctrl, targets, genes)
        row = dict(line=held, method="original_transfer", strength=0., n_targets=len(targets), **scores)
        results.append(row)
        print(row, flush=True)
        out = ROOT / "data/calibration"
        pd.DataFrame(results).to_csv(out / f"{args.name}.csv", index=False)
        (out / f"{args.name}_config.json").write_text(json.dumps(vars(args), indent=2))


if __name__ == "__main__":
    main()
