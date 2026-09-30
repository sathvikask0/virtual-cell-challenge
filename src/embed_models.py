"""Predict knockdown responses from a gene embedding of the knocked-down gene.

Input  x: embedding of the knocked-down gene (e.g. data/embeddings/coexpr.npz)
Output y: share change vs control (18080) and log total ratio vs control (1)
Predicted counts = (control share + share change) * control total * exp(log total ratio)

Models:
  knn    average y of the k training genes with the most similar embedding (cosine)
  ridge  linear regression from x to y with L2 penalty alpha

Hyperparameters (embedding dims, k, alpha) are picked by 5-fold cross-validation over
train genes, then the chosen setting is refit on all train genes and scored on val.
Test is not touched here.

Usage: python src/embed_models.py [embedding name, default coexpr]
"""
import sys
from itertools import product
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

from baseline import load, score

ROOT = Path(__file__).resolve().parent.parent
DIMS = [16, 64, 128]
KS = [1, 3, 5, 10, 20]
ALPHAS = [0.1, 1, 10, 100]


def embed(name, targets, dim):
    d = np.load(ROOT / "data/embeddings" / f"{name}.npz")
    idx = {g: i for i, g in enumerate(d["genes"])}
    X = d["emb"][[idx[g] for g in targets], :dim]
    return X / np.linalg.norm(X, axis=1, keepdims=True)


def fit_predict(model, param, Xtr, Ytr, Xte):
    if model == "mean":
        return np.repeat(Ytr.mean(0, keepdims=True), len(Xte), axis=0)
    if model == "knn":
        nn = np.argsort(-(Xte @ Xtr.T), axis=1)[:, :param]
        return Ytr[nn].mean(1)
    if model == "ridge":
        return Ridge(alpha=param).fit(Xtr, Ytr).predict(Xte)


def to_counts(Y, ctrl):
    share = np.clip(ctrl["share"] + Y[:, :-1], 0, None)
    share /= share.sum(1, keepdims=True)
    return share * (ctrl["total"] * np.exp(Y[:, -1:]))


def targets_y(split):
    targets, share, total, counts, ctrl = load(split)
    Y = np.hstack([share - ctrl["share"], np.log(total / ctrl["total"])[:, None]])
    return targets, Y, counts, ctrl


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "coexpr"
    tr_t, tr_Y, tr_counts, ctrl = targets_y("train")
    settings = ([("mean", None, DIMS[0])]
                + [("knn", k, d) for d, k in product(DIMS, KS)]
                + [("ridge", a, d) for d, a in product(DIMS, ALPHAS)])

    # cross-validation on train genes: every gene gets predicted once while held out
    folds = list(KFold(5, shuffle=True, random_state=0).split(tr_t))
    cv = {}
    for model, param, dim in settings:
        X = embed(name, tr_t, dim)
        pred = np.zeros_like(tr_Y)
        for tr, te in folds:
            pred[te] = fit_predict(model, param, X[tr], tr_Y[tr], X[te])
        cv[(model, param, dim)] = score(to_counts(pred, ctrl), tr_counts, ctrl["counts"])

    print(f"embedding: {name}\n\n5-fold CV on 150 train genes")
    print(f"{'model':6} {'param':>6} {'dims':>5} {'pearson_delta':>14} {'mae':>7} {'total_err':>10} {'discrim':>8}")
    for (model, param, dim), s in cv.items():
        print(f"{model:6} {str(param):>6} {dim:5} {s['pearson_delta']:14.3f} {s['mae']:7.3f} "
              f"{s['total_err']:10.3f} {s['discrimination']:8.3f}")

    # refit best setting per model family on all train genes, score on val
    va_t, _, va_counts, va_ctrl = targets_y("val")
    print("\nval (50 held-out genes), best CV setting per model by pearson_delta")
    print(f"{'model':6} {'param':>6} {'dims':>5} {'pearson_delta':>14} {'mae':>7} {'total_err':>10} {'discrim':>8}")
    for family in ["mean", "knn", "ridge"]:
        model, param, dim = max((k for k in cv if k[0] == family),
                                key=lambda k: cv[k]["pearson_delta"])
        pred = fit_predict(model, param, embed(name, tr_t, dim), tr_Y, embed(name, va_t, dim))
        s = score(to_counts(pred, va_ctrl), va_counts, va_ctrl["counts"])
        print(f"{model:6} {str(param):>6} {dim:5} {s['pearson_delta']:14.3f} {s['mae']:7.3f} "
              f"{s['total_err']:10.3f} {s['discrimination']:8.3f}")


if __name__ == "__main__":
    main()
