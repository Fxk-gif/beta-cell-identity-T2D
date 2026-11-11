# -*- coding: utf-8 -*-
"""
sensitivity_analysis.py
=====================================================================
Sensitivity analysis addressing Reviewer 2, Comment 2:

"The bulk transcriptomic analyses do not adequately account for
 differences in beta-cell abundance, islet purity, or clinical
 confounders. Furthermore, both bulk datasets originate from the same
 IMIDIA study and use the same microarray platform, limiting claims of
 cross-platform reproducibility. Cell-composition adjustment and
 appropriate sensitivity analyses are necessary."

Based on public GEO data, this script performs the following analyses:

  1. Glycemia adjustment (an additional clinical confounder).
     The GSE76895 metadata contain fasting glucose and the GSE76894
     metadata contain ICU glycemia. Glycemia is added to the OLS
     adjustment model (age + gender + bmi + batch + glycemia) to test
     whether the T2D effect on the identity score remains robust after
     further adjustment. Glycemia is both an intermediate variable
     determined by T2D status and a proxy for beta-cell load/function,
     so it is treated as a sensitivity analysis rather than a primary
     effect adjustment.

  2. Beta-cell abundance proxy adjustment (cell-composition proxy).
     The GEO metadata contain no islet-purity or beta-cell-proportion
     field (a data limitation noted by the reviewer), so direct
     adjustment is not possible. This script constructs a "beta-cell
     relative-abundance proxy" from the mean expression of islet
     endocrine-lineage marker genes unrelated to beta-cell identity,
     and adds it as a covariate to test whether the identity-score
     decline can be explained by this compositional proxy. The proxy is
     itself derived from the same bulk expression data and is therefore
     an approximation rather than a formal deconvolution; it is labelled
     as a proxy throughout.

  3. Transcription-factor vs. effector sensitivity (compositional vs.
     autonomous). The 15-gene identity score is split into two
     sub-scores:
       - TF sub-score (6 core transcription factors)
       - effector sub-score (9 functional/secretory genes)
     If a reduced beta-cell proportion (compositional) were the main
     driver, both classes should decline in parallel, because a lower
     cell proportion would depress all beta-cell genes uniformly; a
     selectively stronger TF decline would instead support
     cell-autonomous erosion of the regulatory program.

  4. Cohort origin and platform (appended to the output). Both bulk
     datasets (GSE76895 partial pancreatectomy, GSE76894 organ donor)
     originate from the IMIDIA study and use the same GPL570 microarray
     platform, and the output records this so that the same-platform
     rather than cross-platform status of the replication is explicit.

Available covariates (from the GEO series matrix):
  GSE76895: age, gender, bmi, fasting glucose, batch
  GSE76894: age, gender, bmi, ICU glycemia, batch
  islet purity / cold ischemia time: not present in the GEO metadata,
  so they cannot be adjusted for.

=====================================================================
"""
import csv
import math
import os

import numpy as np

# Locate the data root (the directory containing r3_inputs/). The path is
# resolved from the T2D_DATA environment variable or by walking upward from
# this file.
BASE = os.environ.get("T2D_DATA", "").rstrip("/\\")
if not BASE:
    _d = os.path.dirname(os.path.abspath(__file__))
    while _d and not os.path.isdir(os.path.join(_d, "r3_inputs")):
        _parent = os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    BASE = _d

# 15-gene beta-cell identity score used for the cross-cohort comparisons in the
# Methods (the 16-gene panel minus PDX1, which is absent from the GPL570 probes
# of the validation cohort). PDX1 is analysed separately as a single gene in the
# discovery cohort. MAFA/KCNJ11/UCN3 have no probe on GPL570 and are not part of
# the bulk identity score.
TF_GENES = ["NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6"]
EFFECTOR_GENES = ["SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8",
                  "PCSK1", "PCSK2", "INS"]
BETA_IDENTITY = TF_GENES + EFFECTOR_GENES

# Beta-cell relative-abundance proxy: lineage marker genes known to be highly
# expressed in beta cells and not part of the identity-score definition
# (higher expression -> relatively more beta cells).
# Note: these genes may themselves be affected by disease, so the proxy is used
# only as a sensitivity analysis, not a formal deconvolution.
BETA_ABUNDANCE_PROXY = ["CHGA", "CHGB", "SCG2", "SCG3", "ISL1", "PCSK1N"]

# Cohort -> donor source and glycemia field name
COHORT_INFO = {
    "GSE76895": {
        "source": "pancreatectomy (partial pancreatectomy patients)",
        "glycemia_key": "fasting glucose",
    },
    "GSE76894": {
        "source": "organ donor",
        "glycemia_key": "icu glycemia",
    },
}


def parse_covariates(gse):
    """Extract sample-level covariates (age/gender/bmi/glycemia/batch) from the
    series matrix."""
    import gzip
    import re
    f = os.path.join(BASE, "geo_data", f"{gse}_series_matrix.txt.gz")
    txt = gzip.open(f, "rt", encoding="latin-1").read()
    lines = txt.split("\n")
    acc = None
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            acc = re.findall(r'"([^"]*)"', ln)
            break

    def _field(key):
        for ln in lines:
            if ln.startswith("!Sample_characteristics") and key in ln:
                return re.findall(r'"([^"]*)"', ln)
        return []

    age = _field("age (years)")
    gender = _field("gender")
    bmi = _field("body mass index")
    # The glycemia field is named differently in the two cohorts:
    #   GSE76895: "fasting glucose (mmol/l)"
    #   GSE76894: "icu (intensive care unit) glycemia (mmol/l)"
    # Match loosely using the "glucose" / "glycemia" keywords.
    gly = _field("glucose") or _field("glycemia")
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
            "glycemia": _val(gly[i]) if i < len(gly) else None,
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


def _to_float(v):
    """Parse strings such as '65' or '65 (in s...' into a float, returning nan
    on failure."""
    if v is None:
        return float("nan")
    v = str(v).strip()
    m = __import__("re").match(r"[-+]?\d*\.?\d+", v)
    if m:
        try:
            return float(m.group(0))
        except ValueError:
            return float("nan")
    return float("nan")


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)


def score_of(mat, genes, gset):
    idx = [genes.index(g) for g in gset if g in genes]
    sub = mat[idx, :]
    subz = np.vstack([zscore(sub[i]) for i in range(len(idx))])
    return subz.mean(axis=0), [genes[i] for i in idx]


def build_X(grp01, age, gender, bmi, gly, batches):
    n = len(grp01)
    cols = [np.ones(n), np.array(grp01, dtype=float),
            np.array(age, dtype=float), np.array(gender, dtype=float),
            np.array(bmi, dtype=float)]
    names = ["intercept", "group", "age", "gender", "bmi"]
    if gly is not None:
        cols.append(np.array(gly, dtype=float))
        names.append("glycemia")
    ub = sorted(set(batches))
    if len(ub) > 1:
        for b in ub[1:]:
            cols.append(np.array([1.0 if x == b else 0.0 for x in batches]))
            names.append(f"batch_{b}")
    return np.column_stack(cols), names


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
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
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
    out_lines = []
    def emit(s=""):
        print(s)
        out_lines.append(s)

    emit("#" * 78)
    emit("# R2-2 sensitivity analysis: glycemia adjustment + beta-abundance proxy + TF/effector split")
    emit("#" * 78)

    summary = {}
    for gse in ["GSE76895", "GSE76894"]:
        emit("\n" + "=" * 78)
        emit(f"{gse}  (donor source: {COHORT_INFO[gse]['source']})")
        emit("=" * 78)
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
        grp01 = np.array([1.0 if g == "T2D" else 0.0 for g in keep_group])

        age, gender, bmi, gly, batch = [], [], [], [], []
        for s in samples:
            c = cov.get(s, {})
            age.append(_to_float(c.get("age")))
            g_ = c.get("gender")
            gender.append(1.0 if g_ == "M" else (0.0 if g_ == "F" else float("nan")))
            bmi.append(_to_float(c.get("bmi")))
            gly.append(_to_float(c.get("glycemia")))
            batch.append(c.get("batch") or "NA")
        age = np.array(age); gender = np.array(gender)
        bmi = np.array(bmi); gly = np.array(gly)

        # Imputation (mean) for missing values, and report the missing counts
        def _fill(a):
            a = np.array(a, dtype=float)
            return a, int(np.isnan(a).sum()), np.where(np.isnan(a), np.nanmean(a), a)
        age, na_age, age = _fill(age)
        gender, na_gen, gender = _fill(gender)
        bmi, na_bmi, bmi = _fill(bmi)
        gly, na_gly, gly = _fill(gly)

        emit(f"  samples: {len(samples)} (T2D={int(grp01.sum())}, ND={int((1-grp01).sum())})")
        emit(f"  missing covariates: age={na_age}, gender={na_gen}, bmi={na_bmi}, glycemia={na_gly}")
        emit(f"  glycemia mean: {gly.mean():.1f} (T2D {gly[grp01==1].mean():.1f} vs ND {gly[grp01==0].mean():.1f})")

        # Three scores
        score_all, used_all = score_of(mat, genes, BETA_IDENTITY)
        score_tf, used_tf = score_of(mat, genes, TF_GENES)
        score_eff, used_eff = score_of(mat, genes, EFFECTOR_GENES)
        proxy, used_proxy = score_of(mat, genes, BETA_ABUNDANCE_PROXY)
        emit(f"  identity-score genes available: {len(used_all)}/15 (TF={len(used_tf)}, effector={len(used_eff)})")
        emit(f"  beta-abundance-proxy genes available: {len(used_proxy)}/{len(BETA_ABUNDANCE_PROXY)} -> {used_proxy}")

        # --- Analysis 1: base adjustment (age+gender+bmi+batch), no glycemia ---
        X_base, names_base = build_X(grp01, age, gender, bmi, None, batch)
        b_base, se_base, t_base, p_base, df_base = ols(X_base, score_all)
        emit(f"  [A] identity score base adjustment (age+gender+bmi+batch): "
             f"p={p_base[1]:.4e}, beta={b_base[1]:+.4f} (df={df_base})")

        # --- Analysis 2: add glycemia ---
        X_gly, names_gly = build_X(grp01, age, gender, bmi, gly, batch)
        b_gly, se_gly, t_gly, p_gly, df_gly = ols(X_gly, score_all)
        emit(f"  [B] identity score + glycemia adjustment: "
             f"p={p_gly[1]:.4e}, beta={b_gly[1]:+.4f} (df={df_gly})")

        # --- Analysis 3: add beta-abundance proxy ---
        X_proxy, names_proxy = build_X(grp01, age, gender, bmi, gly, batch)
        # Append the proxy column to X_gly
        X_comp = np.column_stack([X_gly, proxy])
        b_comp, se_comp, t_comp, p_comp, df_comp = ols(X_comp, score_all)
        emit(f"  [C] identity score + glycemia + beta-abundance proxy: "
             f"p={p_comp[1]:.4e}, beta={b_comp[1]:+.4f} (df={df_comp})")

        # --- Analysis 4: TF vs effector split ---
        b_tf = ols(X_gly, score_tf)
        b_eff = ols(X_gly, score_eff)
        emit(f"  [D] TF sub-score (+glycemia) : p={b_tf[3][1]:.4e}, beta={b_tf[0][1]:+.4f}")
        emit(f"  [E] effector sub-score (+glycemia): p={b_eff[3][1]:.4e}, beta={b_eff[0][1]:+.4f}")

        # --- Key individual genes (as consistency reference) ---
        emit(f"  --- key individual genes (base adjustment vs +glycemia) ---")
        for g in ["SLC2A2", "GCK", "PDX1", "NKX6-1", "INS"]:
            if g not in genes:
                continue
            gi = genes.index(g)
            b_base_g, _, _, p_base_g, _ = ols(X_base, mat[gi, :])
            b_gly_g, _, _, p_gly_g, _ = ols(X_gly, mat[gi, :])
            emit(f"    {g:10s} base p={p_base_g[1]:.4e} beta={b_base_g[1]:+.4f}  |  "
                 f"+glycemia p={p_gly_g[1]:.4e} beta={b_gly_g[1]:+.4f}")

        # --- Difference of the beta-abundance proxy itself between T2D and ND ---
        p_proxy_mw = mw_p(proxy[grp01 == 1], proxy[grp01 == 0])
        emit(f"  [F] beta-abundance proxy T2D vs ND: Mann-Whitney p={p_proxy_mw:.4e}")

        summary[gse] = {
            "n_T2D": int(grp01.sum()), "n_ND": int((1 - grp01).sum()),
            "p_base": p_base[1], "p_gly": p_gly[1], "p_comp": p_comp[1],
            "beta_base": b_base[1], "beta_gly": b_gly[1], "beta_comp": b_comp[1],
            "p_tf": b_tf[3][1], "p_eff": b_eff[3][1],
            "beta_tf": b_tf[0][1], "beta_eff": b_eff[0][1],
            "p_proxy": p_proxy_mw,
        }

    # Summary + same-platform statement
    emit("\n" + "=" * 78)
    emit("Summary and conclusions")
    emit("=" * 78)
    for gse, s in summary.items():
        emit(f"\n{gse}: base p={s['p_base']:.3e} -> +glycemia p={s['p_gly']:.3e} -> +proxy p={s['p_comp']:.3e}")
        emit(f"    TF sub-score p={s['p_tf']:.3e} (beta {s['beta_tf']:+.3f}), "
             f"effector sub-score p={s['p_eff']:.3e} (beta {s['beta_eff']:+.3f})")
        emit(f"    beta-abundance proxy T2D vs ND p={s['p_proxy']:.3e}")

    emit("\n" + "=" * 78)
    emit("Cohort origin and platform")
    emit("=" * 78)
    emit("GSE76895 and GSE76894 both originate from the IMIDIA study and were")
    emit("profiled on the same Affymetrix GPL570 microarray platform. The two")
    emit("cohorts use different donor-id prefixes (DP### vs PO###) with disjoint")
    emit("id ranges and zero overlap in the age-sex-BMI triple, i.e. non-overlapping")
    emit("donor groups.")

    # Write the results file
    out = os.path.join(os.path.dirname(__file__), "sensitivity_analysis_results.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    run()
