"""EBM comparator, step 2: fit and predict (docs/ebm_comparator_20261028.md).

Reads the design tables written by tests/ebm_export.R and fits the
Explainable Boosting Machine with its DEFAULT settings (pairwise interactions
included, native missing-value handling) and a fixed seed, exactly as the
`xgb_feat` booster is fitted:

  out of fold   for each fold f, fit on the training rows outside f using
                design_fold<f> (that fold's out-of-fold priors), predict the
                rows in f; seed = SEED + f, as xgb_oof_perfold()
  final         fit on all training rows using design_final (frozen final
                priors); apply to design_eicu and design_test (apply only)

Writes predictions and the final model under --out. Prints aggregates only
(row counts, timings, term counts): the tables are derived from credentialed
data and never leave this machine (CLAUDE.md hard rule 1).

  .venv-ebm/Scripts/python -I python/ebm_fit.py --benchmark
  .venv-ebm/Scripts/python -I python/ebm_fit.py
"""
import argparse
import os
import pickle
import time

import numpy as np
import pandas as pd
from interpret.glassbox import ExplainableBoostingClassifier

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="out/cache/ebm")
ap.add_argument("--out", default="out/cache/ebm/fit")
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--benchmark", action="store_true", help="fit fold 1 only and report its time")
ap.add_argument("--additive", action="store_true",
                help="main effects only (interactions=0): the purely additive jointly fitted model (author, 2026-10-28)")
args = ap.parse_args()
if args.additive and args.out == "out/cache/ebm/fit":
    args.out = "out/cache/ebm/fit_additive"
os.makedirs(args.out, exist_ok=True)

feats = open(os.path.join(args.dir, "feature_names.txt"), encoding="utf-8").read().splitlines()
rows = pd.read_parquet(os.path.join(args.dir, "rows.parquet"))
y = np.asarray(rows["y"])
fold = np.asarray(rows["fold"])


def read(name):
    d = pd.read_parquet(os.path.join(args.dir, name + ".parquet"))
    if list(d.columns) != feats:
        raise SystemExit(f"{name}: column set differs from feature_names.txt")
    return d


def ebm(seed):
    # Defaults throughout; only the seed is set (decision 3). --additive turns
    # the pairwise interactions off and changes nothing else.
    if args.additive:
        return ExplainableBoostingClassifier(random_state=seed, interactions=0)
    return ExplainableBoostingClassifier(random_state=seed)


def report(label, m, t0, n):
    n_int = sum(1 for t in m.term_features_ if len(t) > 1)
    print(f"{label}: {n} training rows, {time.time() - t0:.0f} s, {len(m.term_features_)} terms ({n_int} interactions)", flush=True)


p_oof = np.full(len(y), np.nan)
for f in sorted(set(fold)):
    X = read(f"design_fold{f}")
    if len(X) != len(y):
        raise SystemExit(f"design_fold{f} has {len(X)} rows, expected {len(y)}")
    tr, ho = fold != f, fold == f
    t0 = time.time()
    m = ebm(args.seed + int(f)).fit(X[tr], y[tr])
    report(f"fold {f}", m, t0, int(tr.sum()))
    p_oof[ho] = m.predict_proba(X[ho])[:, 1]
    if args.benchmark:
        raise SystemExit(0)

if np.isnan(p_oof).any():
    raise SystemExit("a training row was never scored out of fold")
pd.DataFrame({"row": rows["row"], "p": p_oof}).to_parquet(os.path.join(args.out, "pred_oof.parquet"))

Xf = read("design_final")
t0 = time.time()
mf = ebm(args.seed).fit(Xf, y)
report("final", mf, t0, len(y))
with open(os.path.join(args.out, "ebm_final.pkl"), "wb") as fh:
    pickle.dump(mf, fh)
for site in ("eicu", "test"):
    Xs = read(f"design_{site}")
    pd.DataFrame({"p": mf.predict_proba(Xs)[:, 1]}).to_parquet(os.path.join(args.out, f"pred_{site}.parquet"))
    print(f"{site}: {len(Xs)} rows predicted", flush=True)
print("done", flush=True)
