#!/usr/bin/env python3
"""scripts/checks/07_bootstrap_ci_crosscheck.py

Non-canonical Python cross-check for scripts/07_bootstrap_ci.R. NO new
modeling: it resamples the already committed out-of-fold predictions. Because
its PR integral and random-number generator are only approximations to the
declared R/PRROC workflow, it writes separate validation files and must never
overwrite the canonical uncertainty tables.

Reproduces, against the manuscript point estimates:
  * AUROC and operating-point metrics exactly,
  * PR-AUC via a Davis-Goadrich integral that matches PRROC::pr.curve(
    ...)$auc.integral to ~0.001-0.002.

Outputs (cross-check only):
  tables/uncertainty/python_crosscheck/bootstrap_ci_python.csv
  tables/uncertainty/python_crosscheck/delta_auroc_prauc_ci_python.csv
  results/uncertainty/delong_tests_python_crosscheck.csv

Method: stratified percentile bootstrap (B=2000, seed=20260620). The leaky
baseline AUROC/PR-AUC point estimate is the mean across the 5 repeated-CV runs;
the bootstrap recomputes that mean within each resample. Operating-point metrics
pool the repeats. Delta AUROC / Delta PR-AUC use a paired stratified bootstrap
(two-sided bootstrap p-value); DeLong is run per repeat (leaky_r vs guarded).
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy import stats

SEED, B, POS, NEG = 20260620, 2000, "pCR", "RD"
rng = np.random.default_rng(SEED)
os.makedirs("tables/uncertainty", exist_ok=True)
os.makedirs("tables/uncertainty/python_crosscheck", exist_ok=True)
os.makedirs("results/uncertainty", exist_ok=True)


def pr_auc_integral(y, score):
    """Davis-Goadrich PR integral (matches PRROC::pr.curve auc.integral)."""
    y = np.asarray(y); score = np.asarray(score, float); P = y.sum()
    if P == 0 or P == len(y):
        return np.nan
    order = np.argsort(-score, kind="mergesort"); s = score[order]; yy = y[order]
    tp = fp = 0; pts = [(0, 0)]; i = 0; n = len(s)
    while i < n:
        j = i
        while j < n and s[j] == s[i]:
            if yy[j] == 1:
                tp += 1
            else:
                fp += 1
            j += 1
        pts.append((tp, fp)); i = j
    area = 0.0
    for (tpa, fpa), (tpb, fpb) in zip(pts[:-1], pts[1:]):
        dtp = tpb - tpa
        if dtp == 0:
            continue
        dfp = fpb - fpa
        for k in range(dtp):
            tp0, fp0 = tpa + k, fpa + dfp * k / dtp
            tp1, fp1 = tpa + k + 1, fpa + dfp * (k + 1) / dtp
            prec0 = tp0 / (tp0 + fp0) if (tp0 + fp0) > 0 else 1.0
            prec1 = tp1 / (tp1 + fp1)
            area += 0.5 * (prec0 + prec1) * ((tp1 - tp0) / P)
    return area


def op_metrics(truth, pred):
    y = (np.asarray(truth) == POS); p = (np.asarray(pred) == POS)
    tp = np.sum(y & p); tn = np.sum(~y & ~p); fp = np.sum(~y & p); fn = np.sum(y & ~p)
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den > 0 else np.nan
    return dict(balanced_accuracy=np.nanmean([sens, spec]), mcc=mcc,
                sensitivity=sens, specificity=spec)


def delong_test(y, s1, s2):
    """Fast DeLong test for two correlated ROC AUCs (Sun & Xu, 2014)."""
    def midrank(x):
        J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
        while i < N:
            j = i
            while j < N and Z[j] == Z[i]:
                j += 1
            T[i:j] = 0.5 * (i + j - 1) + 1
            i = j
        T2 = np.empty(N); T2[J] = T; return T2
    y = np.asarray(y); m = int(y.sum()); n = len(y) - m; pos = y == 1
    preds = np.vstack([s1, s2]); posS = preds[:, pos]; negS = preds[:, ~pos]
    tx = np.array([midrank(posS[r]) for r in range(2)])
    ty = np.array([midrank(negS[r]) for r in range(2)])
    tz = np.array([midrank(np.concatenate([posS[r], negS[r]])) for r in range(2)])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    L = np.array([1, -1]); var = L @ cov @ L
    z = (aucs[0] - aucs[1]) / np.sqrt(var) if var > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if var > 0 else np.nan
    return aucs[0], aucs[1], z, p


# ---- load ----
leaky = pd.read_csv("results/pilot_gse25055/leaky_baseline_predictions.csv")
nested = pd.read_csv("results/pilot_gse25055/nested_smo_svm_predictions.csv")
ext = pd.read_csv("results/external_validation_gse25065/gse25065_external_predictions.csv")
reps = sorted(leaky.repeat_id.unique())
ids = list(nested.sample_id)
truth_map = nested.set_index("sample_id").truth.to_dict()
y_int = np.array([1 if truth_map[i] == POS else 0 for i in ids])
leaky_prob = {r: leaky[leaky.repeat_id == r].set_index("sample_id").prob_pos.to_dict() for r in reps}
leaky_pred = {r: leaky[leaky.repeat_id == r].set_index("sample_id").predicted.to_dict() for r in reps}
nested_prob = nested.set_index("sample_id").prob_pos.to_dict()
nested_pred = nested.set_index("sample_id").predicted.to_dict()
y_ext = (ext.truth.values == POS).astype(int)
METRICS = ["auroc", "pr_auc", "balanced_accuracy", "mcc", "sensitivity", "specificity"]


def auroc(y, s):
    return roc_auc_score(y, s)


def leaky_point():
    a = np.mean([auroc(y_int, np.array([leaky_prob[r][i] for i in ids])) for r in reps])
    pr = np.mean([pr_auc_integral(y_int, np.array([leaky_prob[r][i] for i in ids])) for r in reps])
    return dict(auroc=a, pr_auc=pr, **op_metrics(leaky.truth.values, leaky.predicted.values))


def nested_point():
    s = np.array([nested_prob[i] for i in ids])
    return dict(auroc=auroc(y_int, s), pr_auc=pr_auc_integral(y_int, s),
                **op_metrics([truth_map[i] for i in ids], [nested_pred[i] for i in ids]))


def ext_point():
    return dict(auroc=auroc(y_ext, ext.prob_pos.values),
                pr_auc=pr_auc_integral(y_ext, ext.prob_pos.values),
                **op_metrics(ext.truth.values, ext.predicted.values))


pos_idx = np.where(y_int == 1)[0]; neg_idx = np.where(y_int == 0)[0]
pos_e = np.where(y_ext == 1)[0]; neg_e = np.where(y_ext == 0)[0]


def resample_int():
    return np.concatenate([rng.choice(pos_idx, len(pos_idx), replace=True),
                           rng.choice(neg_idx, len(neg_idx), replace=True)])


def resample_ext():
    return np.concatenate([rng.choice(pos_e, len(pos_e), replace=True),
                           rng.choice(neg_e, len(neg_e), replace=True)])


points = {"leaky": leaky_point(), "nested": nested_point(), "external": ext_point()}
labels = {"leaky": "A_leaky_baseline", "nested": "B_guarded_nested",
          "external": "B_guarded_nested_external"}
cohorts = {"leaky": "GSE25055", "nested": "GSE25055", "external": "GSE25065"}

rows = []
for kind in ["leaky", "nested", "external"]:
    acc = {m: [] for m in METRICS}
    for _ in range(B):
        if kind == "external":
            idx = resample_ext(); yb = y_ext[idx]
            acc["auroc"].append(auroc(yb, ext.prob_pos.values[idx]))
            acc["pr_auc"].append(pr_auc_integral(yb, ext.prob_pos.values[idx]))
            op = op_metrics(ext.truth.values[idx], ext.predicted.values[idx])
        else:
            idx = resample_int(); rids = [ids[k] for k in idx]; yb = y_int[idx]
            if kind == "leaky":
                acc["auroc"].append(np.mean([auroc(yb, np.array([leaky_prob[r][i] for i in rids])) for r in reps]))
                acc["pr_auc"].append(np.mean([pr_auc_integral(yb, np.array([leaky_prob[r][i] for i in rids])) for r in reps]))
                op = op_metrics([truth_map[i] for i in rids] * len(reps),
                                np.concatenate([[leaky_pred[r][i] for i in rids] for r in reps]))
            else:
                s = np.array([nested_prob[i] for i in rids])
                acc["auroc"].append(auroc(yb, s)); acc["pr_auc"].append(pr_auc_integral(yb, s))
                op = op_metrics([truth_map[i] for i in rids], [nested_pred[i] for i in rids])
        for m in METRICS[2:]:
            acc[m].append(op[m])
    for m in METRICS:
        lo, hi = np.nanpercentile(acc[m], [2.5, 97.5])
        rows.append(dict(cohort=cohorts[kind], pipeline=labels[kind], metric=m,
                         point=round(points[kind][m], 4), ci_lo=round(lo, 4),
                         ci_hi=round(hi, 4), n_boot=B,
                         method="stratified percentile bootstrap"))
pd.DataFrame(rows).to_csv(
    "tables/uncertainty/python_crosscheck/bootstrap_ci_python.csv", index=False
)

# paired bootstrap delta
dau, dpr = [], []
for _ in range(B):
    idx = resample_int(); rids = [ids[k] for k in idx]; yb = y_int[idx]
    la = np.mean([auroc(yb, np.array([leaky_prob[r][i] for i in rids])) for r in reps])
    lp = np.mean([pr_auc_integral(yb, np.array([leaky_prob[r][i] for i in rids])) for r in reps])
    dau.append(la - auroc(yb, np.array([nested_prob[i] for i in rids])))
    dpr.append(lp - pr_auc_integral(yb, np.array([nested_prob[i] for i in rids])))
dau, dpr = np.array(dau), np.array(dpr)


def p2(a):
    return round(2 * min(np.mean(a > 0), np.mean(a < 0)), 4)


pd.DataFrame([
    dict(comparison="leaky - guarded (GSE25055)", metric="delta_auroc",
         point=round(points["leaky"]["auroc"] - points["nested"]["auroc"], 4),
         ci_lo=round(np.percentile(dau, 2.5), 4), ci_hi=round(np.percentile(dau, 97.5), 4),
         boot_p=p2(dau), n_boot=B, method="paired stratified bootstrap"),
    dict(comparison="leaky - guarded (GSE25055)", metric="delta_pr_auc",
         point=round(points["leaky"]["pr_auc"] - points["nested"]["pr_auc"], 4),
         ci_lo=round(np.percentile(dpr, 2.5), 4), ci_hi=round(np.percentile(dpr, 97.5), 4),
         boot_p=p2(dpr), n_boot=B, method="paired stratified bootstrap"),
]).to_csv(
    "tables/uncertainty/python_crosscheck/delta_auroc_prauc_ci_python.csv",
    index=False,
)

# DeLong per repeat
drows = []
for r in reps:
    s1 = np.array([leaky_prob[r][i] for i in ids]); s2 = np.array([nested_prob[i] for i in ids])
    a1, a2, z, p = delong_test(y_int, s1, s2)
    drows.append(dict(comparison="leaky_repeat_%d vs guarded_nested" % r,
                      auroc_leaky=round(a1, 4), auroc_guarded=round(a2, 4),
                      auroc_diff=round(a1 - a2, 4), delong_z=round(z, 4),
                      delong_p=round(p, 4), n=len(ids)))
pd.DataFrame(drows).to_csv(
    "results/uncertainty/delong_tests_python_crosscheck.csv", index=False
)
print("[07-crosscheck] wrote separate non-canonical validation files; DeLong median p =",
      round(float(np.median([d["delong_p"] for d in drows])), 4))
