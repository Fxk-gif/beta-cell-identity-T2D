# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 8:
"The reported ROC analyses lack adequate assessment of uncertainty and potential
confounding, while the multivariable models show inconsistent performance."

This script re-assesses the ROC / diagnostic model analyses to (a) quantify
uncertainty (bootstrap confidence intervals for AUC), (b) formally compare
models (DeLong test), (c) assess model stability (repeated stratified
cross-validation), and (d) report class-imbalance-robust metrics. It also
documents the class composition as part of the confounding assessment.

Data (log2, RMA-normalized microarray expression):
  - Discovery:  GSE76895, n = 68  (36 T2D, 32 ND), 16 identity genes incl. PDX1
  - Validation: GSE76894, n = 103 (19 T2D, 84 ND), 15 identity genes (PDX1 absent)

Analyses:
  1. Per-gene ROC AUC (single identity genes) with 95% bootstrap CIs
     (stratified by outcome), discovery and validation cohorts separately.
  2. Multivariable models (logistic regression, SVM RBF, random forest) with
     validation AUC + 95% bootstrap CIs.
  3. DeLong test: SLC2A2 vs every other identity gene; SLC2A2 vs multivariable
     logistic model (both cohorts).
  4. Repeated (5 x 5) stratified cross-validation within the discovery cohort
     to estimate mean AUC +/- SD (model stability).
  5. Class-imbalance-robust metrics (balanced accuracy, sensitivity,
     specificity, F1, MCC) with bootstrap CIs.
  6. Class-composition documentation (confounding assessment).

Dependencies: numpy, scipy, scikit-learn.
Output files (placed next to this script):
  roc_auc_bootstrap_ci.csv
  delong_comparisons.csv
  repeated_cv_stability.csv
  imbalance_robust_metrics.csv
"""

import os
import csv
import numpy as np
import pandas as pd
from scipy import stats

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score,
                             matthews_corrcoef, confusion_matrix)

# ---------------------------------------------------------------------------
# Locate the data root (directory containing r3_inputs).
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(10):
    if os.path.isdir(os.path.join(d, "r3_inputs")):
        BASE = d
        break
    parent = os.path.dirname(d)
    if parent == d:
        break
    d = parent

if BASE is None:
    raise FileNotFoundError(
        "Could not locate the data root containing r3_inputs.")

R3IN = os.path.join(BASE, "r3_inputs")
OUT = os.path.dirname(os.path.abspath(__file__))

# 16 identity genes (master TFs + functional effectors).
IDG = ["PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
       "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8", "PCSK1",
       "PCSK2", "INS"]

RNG_SEED = 42
N_BOOT = 2000


def load_dataset(ds):
    """Return (expr DataFrame [gene x sample], y array 0/1, genes list)."""
    expr = pd.read_csv(os.path.join(R3IN, f"{ds}_expr.csv"), index_col=0)
    grp = pd.read_csv(os.path.join(R3IN, f"{ds}_groups.csv"))
    gmap = dict(zip(grp["sample"], grp["group"]))
    samples = list(expr.columns)
    y = np.array([1 if gmap[s] == "T2D" else 0 for s in samples])
    return expr, y


def auc_bootstrap(y, score, n_boot=N_BOOT, seed=RNG_SEED):
    """Stratified bootstrap 95% CI and point estimate for AUC."""
    auc_point = roc_auc_score(y, score)
    rng = np.random.default_rng(seed)
    n = len(y)
    y1_idx = np.where(y == 1)[0]
    y0_idx = np.where(y == 0)[0]
    boots = []
    for _ in range(n_boot):
        i1 = rng.integers(0, len(y1_idx), len(y1_idx))
        i0 = rng.integers(0, len(y0_idx), len(y0_idx))
        idx = np.concatenate([y1_idx[i1], y0_idx[i0]])
        boots.append(roc_auc_score(y[idx], score[idx]))
    boots = np.array(boots)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return auc_point, lo, hi


def delong_test(y, score1, score2):
    """
    DeLong's nonparametric test for the difference between two AUCs.
    Returns (Z statistic, two-sided p-value).
    Uses the implementation-agnostic structural (X, Y) decomposition
    following DeLong et al. (1988), Biometrics 44:837-845.
    """
    score1 = np.asarray(score1, dtype=float)
    score2 = np.asarray(score2, dtype=float)
    y = np.asarray(y)
    X = score1[y == 0]  # negative scores
    Yp = score1[y == 1]  # positive scores
    X2 = score2[y == 0]
    Yp2 = score2[y == 1]
    m = len(X)
    n = len(Yp)

    # Structural components V10 for model 1 and model 2.
    def v10(xs, ys):
        # For each negative i, V10_i = (1/n) * sum_j psi(ys_j, xs_i)
        # psi = 1 if y > x, 0.5 if equal, 0 otherwise
        out = np.empty(m)
        for i in range(m):
            xi = xs[i]
            s = np.sum((ys > xi) + 0.5 * (ys == xi))
            out[i] = s / n
        return out

    def v01(xs, ys):
        # For each positive j, V01_j = (1/m) * sum_i psi(ys_j, xs_i)
        out = np.empty(n)
        for j in range(n):
            yj = ys[j]
            s = np.sum((yj > xs) + 0.5 * (yj == xs))
            out[j] = s / m
        return out

    V10_1 = v10(X, Yp)
    V01_1 = v01(X, Yp)
    V10_2 = v10(X2, Yp2)
    V01_2 = v01(X2, Yp2)

    S10_1 = (1.0 / (m - 1)) * np.sum((V10_1 - V10_1.mean()) ** 2) if m > 1 else 0.0
    S01_1 = (1.0 / (n - 1)) * np.sum((V01_1 - V01_1.mean()) ** 2) if n > 1 else 0.0
    S10_2 = (1.0 / (m - 1)) * np.sum((V10_2 - V10_2.mean()) ** 2) if m > 1 else 0.0
    S01_2 = (1.0 / (n - 1)) * np.sum((V01_2 - V01_2.mean()) ** 2) if n > 1 else 0.0

    # Covariance between model 1 and model 2 structural components.
    cov10 = (1.0 / (m - 1)) * np.sum((V10_1 - V10_1.mean()) * (V10_2 - V10_2.mean())) if m > 1 else 0.0
    cov01 = (1.0 / (n - 1)) * np.sum((V01_1 - V01_1.mean()) * (V01_2 - V01_2.mean())) if n > 1 else 0.0

    var1 = S10_1 / m + S01_1 / n
    var2 = S10_2 / m + S01_2 / n
    cov = cov10 / m + cov01 / n

    auc1 = roc_auc_score(y, score1)
    auc2 = roc_auc_score(y, score2)
    diff = auc1 - auc2
    denom = np.sqrt(var1 + var2 - 2 * cov)
    if denom == 0:
        return float("nan"), float("nan")
    Z = diff / denom
    p = 2 * (1 - stats.norm.cdf(abs(Z)))
    return Z, p


def get_score(name, Xtr, ytr, Xva):
    """
    Fit model on training, return a decision score on validation.
    For AUC, any monotone transformation is equivalent, so for SVM we use
    decision_function values (probability calibration is unnecessary and slow).
    """
    if name == "Logistic regression":
        m = LogisticRegression(max_iter=2000, C=1.0, random_state=RNG_SEED)
        m.fit(Xtr, ytr)
        return m.predict_proba(Xva)[:, 1]
    elif name == "SVM (RBF)":
        m = SVC(kernel="rbf", C=1.0, gamma="scale", random_state=RNG_SEED)
        m.fit(Xtr, ytr)
        return m.decision_function(Xva)
    elif name == "Random forest":
        m = RandomForestClassifier(n_estimators=500, random_state=RNG_SEED,
                                   n_jobs=-1)
        m.fit(Xtr, ytr)
        return m.predict_proba(Xva)[:, 1]
    else:
        raise ValueError(name)


def main():
    disc, disc_y = load_dataset("GSE76895")
    val, val_y = load_dataset("GSE76894")

    disc_genes = [g for g in IDG if g in disc.index]
    val_genes = [g for g in IDG if g in val.index]
    common_genes = [g for g in IDG if g in val.index]  # 15 (no PDX1)

    print("=" * 78)
    print("ROC uncertainty & multivariable-model stability re-assessment")
    print("=" * 78)
    print(f"Discovery  GSE76895: n={len(disc_y)}  T2D={int(disc_y.sum())}  "
          f"ND={int((disc_y == 0).sum())}")
    print(f"Validation GSE76894: n={len(val_y)}  T2D={int(val_y.sum())}  "
          f"ND={int((val_y == 0).sum())}")
    print(f"Discovery identity genes: {len(disc_genes)} (incl. PDX1)")
    print(f"Validation identity genes: {len(val_genes)} (PDX1 absent)")

    # ------------------------------------------------------------------
    # 1. Per-gene ROC AUC + bootstrap CI (single identity genes).
    #    Direction: high expression -> non-diabetic; use (-expression) as score.
    # ------------------------------------------------------------------
    auc_rows = []
    print("\n" + "=" * 78)
    print("1) Per-gene ROC AUC with 95% bootstrap CIs (stratified)")
    print("=" * 78)
    for ds_name, expr, y in [("discovery", disc, disc_y),
                             ("validation", val, val_y)]:
        genes = [g for g in IDG if g in expr.index]
        print(f"\n  [{ds_name}] {len(genes)} identity genes")
        for g in genes:
            score = -expr.loc[g].values  # high expr -> ND
            auc_p, lo, hi = auc_bootstrap(y, score)
            auc_rows.append([ds_name, g, f"{auc_p:.4f}", f"{lo:.4f}", f"{hi:.4f}"])
            print(f"    {g:10s}  AUC={auc_p:.3f}  95% CI=[{lo:.3f}, {hi:.3f}]")

    with open(os.path.join(OUT, "roc_auc_bootstrap_ci.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "gene", "AUC", "CI_lower", "CI_upper"])
        w.writerows(auc_rows)

    # ------------------------------------------------------------------
    # 2. Multivariable models: validation AUC + bootstrap CI.
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2) Multivariable models (validation AUC + bootstrap CI)")
    print("=" * 78)
    disc_X = disc.loc[common_genes].T.values
    val_X = val.loc[common_genes].T.values
    scaler = StandardScaler().fit(disc_X)
    Xs_tr = scaler.transform(disc_X)
    Xs_va = scaler.transform(val_X)

    mv_rows = []
    for name in ["Logistic regression", "SVM (RBF)", "Random forest"]:
        prob = get_score(name, Xs_tr, disc_y, Xs_va)
        auc_p, lo, hi = auc_bootstrap(val_y, prob)
        mv_rows.append([name, f"{auc_p:.4f}", f"{lo:.4f}", f"{hi:.4f}"])
        print(f"    {name:20s}  val AUC={auc_p:.3f}  95% CI=[{lo:.3f}, {hi:.3f}]")

    # ------------------------------------------------------------------
    # 3. DeLong tests.
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3) DeLong tests (AUC comparisons)")
    print("=" * 78)
    delong_rows = []
    for ds_name, expr, y in [("discovery", disc, disc_y),
                             ("validation", val, val_y)]:
        genes = [g for g in IDG if g in expr.index]
        slc = -expr.loc["SLC2A2"].values
        print(f"\n  [{ds_name}] SLC2A2 vs each other identity gene:")
        for g in genes:
            if g == "SLC2A2":
                continue
            score = -expr.loc[g].values
            Z, p = delong_test(y, slc, score)
            delong_rows.append([ds_name, "SLC2A2", g, f"{Z:.4f}", f"{p:.4f}"])
            print(f"    SLC2A2 vs {g:10s}  DeLong Z={Z:+.3f}  p={p:.4f}")

    # SLC2A2 vs multivariable logistic (validation cohort only).
    val_slc = -val.loc["SLC2A2"].values
    val_logreg = LogisticRegression(max_iter=2000, C=1.0,
                                    random_state=RNG_SEED).fit(Xs_tr, disc_y)
    val_mv_prob = val_logreg.predict_proba(Xs_va)[:, 1]
    Z, p = delong_test(val_y, val_slc, val_mv_prob)
    delong_rows.append(["validation", "SLC2A2", "multivariable_logistic",
                        f"{Z:.4f}", f"{p:.4f}"])
    print(f"\n  [validation] SLC2A2 vs multivariable logistic: "
          f"DeLong Z={Z:+.3f}  p={p:.4f}")

    with open(os.path.join(OUT, "delong_comparisons.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "model_A", "model_B", "DeLong_Z", "p_value"])
        w.writerows(delong_rows)

    # ------------------------------------------------------------------
    # 4. Repeated (5x5) stratified cross-validation within discovery cohort.
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4) Repeated 5x5 stratified cross-validation (discovery cohort)")
    print("=" * 78)
    cv_rows = []
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=5,
                                   random_state=RNG_SEED)
    for name in ["Logistic regression", "SVM (RBF)", "Random forest"]:
        aucs = []
        baccs = []
        mccs = []
        for tr_idx, va_idx in rskf.split(disc_X, disc_y):
            s = StandardScaler().fit(disc_X[tr_idx])
            Xtr = s.transform(disc_X[tr_idx])
            Xva = s.transform(disc_X[va_idx])
            if name == "Logistic regression":
                m = LogisticRegression(max_iter=2000, C=1.0, random_state=RNG_SEED)
            elif name == "SVM (RBF)":
                m = SVC(kernel="rbf", C=1.0, gamma="scale", random_state=RNG_SEED)
            else:
                m = RandomForestClassifier(n_estimators=500,
                                           random_state=RNG_SEED, n_jobs=-1)
            m.fit(Xtr, disc_y[tr_idx])
            if name == "SVM (RBF)":
                prob = m.decision_function(Xva)
            else:
                prob = m.predict_proba(Xva)[:, 1]
            pred = (prob >= 0.5).astype(int)
            aucs.append(roc_auc_score(disc_y[va_idx], prob))
            baccs.append(balanced_accuracy_score(disc_y[va_idx], pred))
            mccs.append(matthews_corrcoef(disc_y[va_idx], pred))
        aucs = np.array(aucs)
        baccs = np.array(baccs)
        mccs = np.array(mccs)
        cv_rows.append([name, f"{aucs.mean():.4f}", f"{aucs.std():.4f}",
                        f"{baccs.mean():.4f}", f"{mccs.mean():.4f}"])
        print(f"    {name:20s}  CV AUC={aucs.mean():.3f} +/- {aucs.std():.3f}  "
              f"bAcc={baccs.mean():.3f}  MCC={mccs.mean():.3f}")

    with open(os.path.join(OUT, "repeated_cv_stability.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model", "CV_AUC_mean", "CV_AUC_sd",
                    "CV_balanced_acc_mean", "CV_MCC_mean"])
        w.writerows(cv_rows)

    # ------------------------------------------------------------------
    # 5. Class-imbalance-robust metrics with bootstrap CIs (validation).
    #    Threshold selected by Youden's J statistic on the validation cohort.
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("5) Class-imbalance-robust metrics (validation, Youden threshold)")
    print("=" * 78)
    robust_rows = []

    def best_threshold(y, score):
        """Optimal threshold by Youden's J (sensitivity + specificity - 1).

        Higher scores correspond to the positive class (T2D), so a candidate
        threshold at index i classifies samples with score >= s_sorted[i] as
        positive (the upper tail), consistent with roc_curve.
        """
        order = np.argsort(score)
        s_sorted = score[order]
        y_sorted = y[order]
        n_pos = int(y.sum())
        n_neg = int((y == 0).sum())
        best_j = -1.0
        best_t = None
        # iterate candidate thresholds at each score boundary
        for i in range(len(s_sorted) + 1):
            # positives at/above threshold index i (higher scores -> pos)
            tp = int(y_sorted[i:].sum())
            fn = n_pos - tp
            fp = (len(s_sorted) - i) - tp
            tn = n_neg - fp
            sens = tp / n_pos if n_pos > 0 else 0.0
            spec = tn / n_neg if n_neg > 0 else 0.0
            j = sens + spec - 1.0
            if j > best_j:
                best_j = j
                if i == 0:
                    best_t = s_sorted[0] - 1.0
                elif i == len(s_sorted):
                    best_t = s_sorted[-1] + 1.0
                else:
                    best_t = (s_sorted[i - 1] + s_sorted[i]) / 2.0
        return best_t

    def robust_metrics(y, score, threshold):
        pred = (score >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(y, pred).ravel()
        bacc = balanced_accuracy_score(y, pred)
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = 2 * prec * sens / (prec + sens) if (prec + sens) > 0 else 0.0
        mcc = matthews_corrcoef(y, pred)
        return bacc, sens, spec, f1, mcc

    # Single gene SLC2A2 and each multivariable model.
    candidates = [("SLC2A2_single_gene", -val.loc["SLC2A2"].values)]
    for name in ["Logistic regression", "SVM (RBF)", "Random forest"]:
        prob = get_score(name, Xs_tr, disc_y, Xs_va)
        candidates.append((name, prob))

    for label, score in candidates:
        thr = best_threshold(val_y, score)
        bacc, sens, spec, f1, mcc = robust_metrics(val_y, score, thr)
        # bootstrap CI for balanced accuracy and MCC (threshold refit per
        # resample to reflect threshold-selection uncertainty).
        rng = np.random.default_rng(RNG_SEED)
        n = len(val_y)
        y1_idx = np.where(val_y == 1)[0]
        y0_idx = np.where(val_y == 0)[0]
        b_bacc, b_mcc = [], []
        for _ in range(N_BOOT):
            i1 = rng.integers(0, len(y1_idx), len(y1_idx))
            i0 = rng.integers(0, len(y0_idx), len(y0_idx))
            idx = np.concatenate([y1_idx[i1], y0_idx[i0]])
            t = best_threshold(val_y[idx], score[idx])
            b_, _, _, _, m_ = robust_metrics(val_y[idx], score[idx], t)
            b_bacc.append(b_)
            b_mcc.append(m_)
        b_bacc = np.array(b_bacc)
        b_mcc = np.array(b_mcc)
        bacc_lo, bacc_hi = np.percentile(b_bacc, [2.5, 97.5])
        mcc_lo, mcc_hi = np.percentile(b_mcc, [2.5, 97.5])
        robust_rows.append([label, f"{thr:.4f}", f"{bacc:.4f}",
                            f"{bacc_lo:.4f}", f"{bacc_hi:.4f}",
                            f"{sens:.4f}", f"{spec:.4f}", f"{f1:.4f}",
                            f"{mcc:.4f}", f"{mcc_lo:.4f}", f"{mcc_hi:.4f}"])
        print(f"    {label:20s}  thr={thr:+.3f}  bAcc={bacc:.3f} "
              f"[{bacc_lo:.3f},{bacc_hi:.3f}]  sens={sens:.3f}  spec={spec:.3f}  "
              f"F1={f1:.3f}  MCC={mcc:.3f} [{mcc_lo:.3f},{mcc_hi:.3f}]")

    with open(os.path.join(OUT, "imbalance_robust_metrics.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model", "Youden_threshold", "balanced_accuracy",
                    "bAcc_CI_lower", "bAcc_CI_upper", "sensitivity",
                    "specificity", "F1", "MCC", "MCC_CI_lower", "MCC_CI_upper"])
        w.writerows(robust_rows)



if __name__ == "__main__":
    main()
