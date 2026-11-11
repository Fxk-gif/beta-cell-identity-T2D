# -*- coding: utf-8 -*-
"""
confounder_adjustment.py
=====================================================================
Confounder-adjustment analysis.

This script extracts sample-level covariates from the GEO series matrices and:
  1. Fits a multiple linear regression for the beta-cell identity score:
        score ~ group(T2D/ND) + age + gender + bmi + batch
     and reports the adjusted coefficient / SE / t / p for the group term
     (OLS normal equations).
  2. Applies the same adjustment to the key single genes and reports p-values
     before and after adjustment.
  3. For the confounder of organ-donor vs pancreatectomy status: GSE76895 is a
     partial-pancreatectomy cohort and GSE76894 is an organ-donor cohort. They
     are different cohorts, so donor source is inherently confounded with cohort
     and cannot be separated as a covariate within a single cohort. This script
     therefore adjusts within each cohort separately and notes this as an
     inherent limitation of the cross-cohort comparison.

Covariates available (from GEO series matrix !Sample_characteristics):
  GSE76895 (PPP partial pancreatectomy): age, gender, bmi, batch
  GSE76894 (OD organ donor):             age, gender, bmi, batch
Cold ischaemia time and islet purity are absent from the GEO metadata and
cannot be adjusted for; this is stated in the response letter.

All results are computed from publicly available data.
=====================================================================
"""
import gzip
import re
import csv
import math
import os

import numpy as np

# Data root (overridable via environment variable). Otherwise walk upward to
# locate the project root containing r3_inputs/.
def _find_data_root():
    env = os.environ.get("T2D_DATA", "").rstrip("/\\")
    if env and os.path.isdir(env):
        return env
    here = os.path.abspath(os.path.dirname(__file__))
    d = here
    for _ in range(12):
        if os.path.isdir(os.path.join(d, "r3_inputs")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return here

BASE = _find_data_root()

# The 15-gene beta-cell identity score (bulk level, uniform panel across
# cohorts). The original 16-gene panel included PDX1, but PDX1 has no probe on
# the GPL570 platform of the validation cohort GSE76894 (208559_at /
# 210937_s_at / 210938_at are all absent). To keep the two cohorts on an
# identical gene panel, the cross-cohort identity score uses 15 genes
# (16 minus PDX1). PDX1 is analyzed separately as a single gene in the
# discovery cohort (retained in KEY_GENES below).
BETA_IDENTITY = ["NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
                 "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8",
                 "PCSK1", "PCSK2", "INS"]
KEY_GENES = ["INS", "SLC2A2", "GCK", "PDX1", "NKX6-1", "MAFB", "IAPP", "SLC30A8"]

# Cohort -> donor-source description
COHORT_SOURCE = {
    "GSE76895": "pancreatectomy (partial pancreatectomy patients)",
    "GSE76894": "organ donor",
}


def parse_covariates(gse):
    f = os.path.join(BASE, "geo_data", f"{gse}_series_matrix.txt.gz")
    txt = gzip.open(f, "rt", encoding="latin-1").read()
    lines = txt.split("\n")
    acc = None
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            acc = [x for x in re.findall(r'"([^"]*)"', ln)]
            break

    def _field(key):
        for ln in lines:
            if ln.startswith("!Sample_characteristics") and key in ln:
                return [x for x in re.findall(r'"([^"]*)"', ln)]
        return []

    age = _field("age (years)")
    gender = _field("gender")
    bmi = _field("body mass index")
    batch = _field("batch")

    def _val(v):
        if v is None:
            return None
        m = re.search(r":\s*(.+)$", v)
        return m.group(1).strip() if m else v.strip()

    out = {}
    for i, g in enumerate(acc):
        out[g] = {
            "age": _val(age[i]) if i < len(age) else None,
            "gender": _val(gender[i]) if i < len(gender) else None,
            "bmi": _val(bmi[i]) if i < len(bmi) else None,
            "batch": _val(batch[i]) if i < len(batch) else None,
        }
    return out


def load_expr(gse):
    path = os.path.join(BASE, "r3_inputs", f"{gse}_expr.csv")
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    return genes, header[1:], mat


def load_groups(gse):
    path = os.path.join(BASE, "r3_inputs", f"{gse}_groups.csv")
    gmap = {}
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gmap[row["sample"]] = row["group"]
    return gmap


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)


def identity_score(mat, genes, gset):
    idx = [genes.index(g) for g in gset if g in genes]
    sub = mat[idx, :]
    subz = np.vstack([zscore(sub[i]) for i in range(len(idx))])
    return subz.mean(axis=0), [genes[i] for i in idx]


def build_X(group, age, gender, bmi, batches):
    n = len(group)
    cols = [np.ones(n), np.array(group, dtype=float),
            np.array(age, dtype=float), np.array(gender, dtype=float),
            np.array(bmi, dtype=float)]
    ub = sorted(set(batches))
    if len(ub) > 1:
        for b in ub[1:]:
            cols.append(np.array([1.0 if x == b else 0.0 for x in batches]))
    return np.column_stack(cols), ub


def _betainc(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lnB = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    def f(t):
        if t <= 0 or t >= 1:
            return 0.0
        return math.exp((a - 1) * math.log(t) + (b - 1) * math.log(1 - t) - lnB)
    N = 2000
    h = x / N
    s = f(0.0) + f(x)
    for i in range(1, N):
        s += (4 if i % 2 == 1 else 2) * f(i * h)
    return s * h / 3.0


def _t_cdf(t, df):
    x = df / (df + t * t)
    return 1 - 0.5 * _betainc(df / 2.0, 0.5, x)


def ols(X, y):
    n, k = X.shape
    XtX = X.T @ X
    Xty = X.T @ y
    beta = np.linalg.solve(XtX, Xty)
    resid = y - X @ beta
    df = n - k
    sigma2 = (resid @ resid) / df
    cov = sigma2 * np.linalg.inv(XtX)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = np.array([2 * (1 - _t_cdf(abs(tv), df)) for tv in t])
    return beta, se, t, p, df


def mw_p(a, b):
    """Two-sided Mann-Whitney U p-value, using the exact null distribution when
    available, and falling back to the asymptotic normal approximation only if
    scipy is unavailable."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    try:
        from scipy import stats as _stats
        # exact method (ties handled)
        return float(_stats.mannwhitneyu(a, b, alternative="two-sided",
                                         method="exact").pvalue)
    except Exception:
        pass
    na, nb = len(a), len(b)
    allv = np.concatenate([a, b])
    ranks = np.argsort(np.argsort(allv)) + 1.0
    ra = ranks[:na]
    U = ra.sum() - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    _, cnt = np.unique(allv, return_counts=True)
    N = na + nb
    sigma2 = na * nb / 12.0 * ((N + 1) - sum((c**3 - c) for c in cnt) / (N * (N - 1)))
    sigma2 = max(sigma2, 1e-12)
    z = (U - mu) / math.sqrt(sigma2)
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def run():
    print("#" * 78)
    print("# R1-5 confounder adjustment: identity score ~ group + age + gender + bmi + batch")
    print("#" * 78)
    summary = {}
    for gse in ["GSE76895", "GSE76894"]:
        print("\n" + "=" * 78)
        print(gse, f"  (donor source: {COHORT_SOURCE[gse]})")
        print("=" * 78)
        genes, samples, mat = load_expr(gse)
        gmap = load_groups(gse)
        cov = parse_covariates(gse)
        keep_idx, keep_group = [], []
        for j, s in enumerate(samples):
            if s in gmap and gmap[s] in ("T2D", "ND"):
                keep_idx.append(j)
                keep_group.append(gmap[s])
        mat = mat[:, keep_idx]
        samples = [samples[j] for j in keep_idx]
        group = keep_group
        age, gender, bmi, batch = [], [], [], []
        for s in samples:
            c = cov.get(s, {})
            a, g_, b, ba = c.get("age"), c.get("gender"), c.get("bmi"), c.get("batch")
            age.append(float(a) if a and a.replace(".", "").replace("-", "").isdigit() else float("nan"))
            gender.append(1.0 if g_ == "M" else (0.0 if g_ == "F" else float("nan")))
            bmi.append(float(b) if b and b.replace(".", "").isdigit() else float("nan"))
            batch.append(ba if ba else "NA")
        age = np.array(age); gender = np.array(gender); bmi = np.array(bmi)
        n_missing = int(np.isnan(age).sum() + np.isnan(gender).sum() + np.isnan(bmi).sum())
        age = np.where(np.isnan(age), np.nanmean(age), age)
        gender = np.where(np.isnan(gender), np.nanmean(gender), gender)
        bmi = np.where(np.isnan(bmi), np.nanmean(bmi), bmi)
        grp01 = np.array([1.0 if g == "T2D" else 0.0 for g in group])
        print(f"  samples: {len(samples)} (T2D={int(grp01.sum())}, ND={int((1-grp01).sum())})")
        print(f"  missing covariate values imputed: {n_missing}")
        print(f"  mean age: {age.mean():.1f} (T2D {age[grp01==1].mean():.1f} vs ND {age[grp01==0].mean():.1f})")
        print(f"  mean BMI: {bmi.mean():.1f} (T2D {bmi[grp01==1].mean():.1f} vs ND {bmi[grp01==0].mean():.1f})")
        print(f"  gender (male proportion): {gender.mean():.2f}")
        score, used = identity_score(mat, genes, BETA_IDENTITY)
        p_unadj = mw_p(score[grp01 == 1], score[grp01 == 0])
        X, ub = build_X(grp01, age, gender, bmi, batch)
        beta, se, t, p_adj, df = ols(X, score)
        print(f"  [identity score] unadjusted Mann-Whitney p = {p_unadj:.4e}")
        print(f"  [identity score] adjusted OLS: group beta={beta[1]:+.4f}, SE={se[1]:.4f}, t={t[1]:+.3f}, p={p_adj[1]:.4e} (df={df})")
        summary[gse] = {"n_T2D": int(grp01.sum()), "n_ND": int((1-grp01).sum()),
                        "score_unadj_p": p_unadj, "score_adj_p": p_adj[1], "score_beta": beta[1]}
        print(f"  --- key genes, before vs after adjustment ---")
        print(f"  {'gene':10s} {'unadj p':>12s} {'adj p':>12s} {'adj beta':>12s}")
        gene_rows = []
        for g in KEY_GENES:
            if g not in genes:
                continue
            gi = genes.index(g)
            gv = mat[gi, :]
            p_u = mw_p(gv[grp01 == 1], gv[grp01 == 0])
            b2, s2, t2, p2, _ = ols(X, gv)
            print(f"  {g:10s} {p_u:12.4e} {p2[1]:12.4e} {b2[1]:+12.4f}")
            gene_rows.append((g, p_u, p2[1], b2[1]))
        summary[gse]["genes"] = gene_rows

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "confounder_adjustment_results.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write("# R1-5 confounder adjustment results (identity score ~ group + age + gender + bmi + batch)\n")
        f.write("# Method: OLS multiple linear regression; the group coefficient is the T2D vs ND difference\n")
        f.write("# Covariate source: GEO series matrix (!Sample_characteristics): age, gender, bmi, batch\n")
        f.write("# Cold ischaemia time / islet purity: absent from GEO metadata, cannot be adjusted (stated in the response letter)\n\n")
        for gse, s in summary.items():
            f.write(f"## {gse} (donor source: {COHORT_SOURCE[gse]})\n")
            f.write(f"  samples: T2D={s['n_T2D']}, ND={s['n_ND']}\n")
            f.write(f"  identity score: unadjusted p={s['score_unadj_p']:.4e}, adjusted p={s['score_adj_p']:.4e}, beta={s['score_beta']:+.4f}\n")
            f.write(f"  key genes (unadjusted p, adjusted p, beta):\n")
            for g, pu, pa, bt in s["genes"]:
                f.write(f"    {g}: {pu:.4e} -> {pa:.4e} (beta {bt:+.4f})\n")
            f.write("\n")


if __name__ == "__main__":
    run()
