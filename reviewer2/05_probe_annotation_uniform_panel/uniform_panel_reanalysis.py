# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 5:
Clarify PDX1 probe annotations, reconcile Figure 3 vs Supplementary Table S5,
and repeat the identity-score analysis using an identical gene panel across
cohorts.

The reviewer correctly noted that (i) PDX1 is unavailable in the validation
cohort (GSE76894), and (ii) PDX1 coverage was inconsistent between Figure 3 and
Supplementary Table S5. We therefore:

  1. Document the GPL570 PDX1 probe annotation and its presence/absence in each
     cohort (read directly from the raw series-matrix files).
  2. Recompute the beta-cell identity score using a uniform 15-gene panel
     (the 16-gene panel minus PDX1) in BOTH cohorts, so the panel is identical
     across cohorts.
  3. Quantify the impact of dropping PDX1 on the discovery-cohort score.

Identity score definition:
    z-score each gene across samples within a cohort, then average the z-scores
    of the panel genes per sample. Group comparison uses a two-sided exact
    Mann-Whitney U test.
"""

import os
import csv
import gzip
import re
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Locate the data root (directory containing r3_inputs and geo_data).
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(10):
    if os.path.isdir(os.path.join(d, "r3_inputs")) and os.path.isdir(os.path.join(d, "geo_data")):
        BASE = d
        break
    parent = os.path.dirname(d)
    if parent == d:
        break
    d = parent

if BASE is None:
    BASE = os.environ.get("T2D_DATA", "")
    if not BASE or not os.path.isdir(os.path.join(BASE, "r3_inputs")):
        raise FileNotFoundError(
            "Could not locate the data root. Set T2D_DATA to the project directory."
        )

R3IN = os.path.join(BASE, "r3_inputs")
GEOD = os.path.join(BASE, "geo_data")

# 16-gene identity panel. PDX1 is removed to form the uniform 15-gene panel
# used across cohorts.
IDG16 = [
    "PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
    "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8", "PCSK1",
    "PCSK2", "INS",
]
IDG15 = [g for g in IDG16 if g != "PDX1"]

# GPL570 (HG-U133 Plus 2.0) PDX1 probes.
PDX1_PROBES = ["208559_at", "210937_s_at", "210938_at"]


def load_expr(ds):
    genes = []
    with open(os.path.join(R3IN, f"{ds}_expr.csv"), encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    gmap = {}
    with open(os.path.join(R3IN, f"{ds}_groups.csv"), encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gmap[row["sample"]] = row["group"]
    return genes, header[1:], mat, gmap


def probe_presence(ds):
    """Check PDX1 probe presence in the raw series-matrix file."""
    path = os.path.join(GEOD, f"{ds}_series_matrix.txt.gz")
    txt = gzip.open(path, "rt", encoding="latin-1").read()
    return {p: (f'"{p}"' in txt) for p in PDX1_PROBES}


def identity_score(ds, panel):
    """z-score average identity score for the given panel in one cohort.

    Returns (scores, groups) aligned to samples kept to T2D/ND.
    """
    genes, samples, mat, gmap = load_expr(ds)
    keep_idx = [j for j, s in enumerate(samples) if gmap.get(s) in ("T2D", "ND")]
    groups = np.array([gmap[samples[j]] for j in keep_idx])
    m = mat[:, keep_idx]

    avail = [g for g in panel if g in genes]
    missing = [g for g in panel if g not in genes]
    gi = [genes.index(g) for g in avail]
    sub = m[gi].astype(float)

    # z-score each gene across the kept samples (ddof=1), then average.
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True, ddof=1)
    z = np.zeros_like(sub, dtype=float)
    ok = sd[:, 0] > 0
    z[ok] = (sub[ok] - mu[ok]) / sd[ok]
    score = z.mean(axis=0)
    return score, groups, avail, missing


def report(ds, panel, label):
    score, groups, avail, missing = identity_score(ds, panel)
    t2d = score[groups == "T2D"]
    nd = score[groups == "ND"]
    # Exact Mann-Whitney U (two-sided), not the asymptotic normal approximation.
    u, p = stats.mannwhitneyu(t2d, nd, alternative="two-sided", method="exact")
    print(f"  [{label}] {ds}")
    print(f"    genes used: {len(avail)}  missing: {missing if missing else 'none'}")
    print(f"    T2D (n={len(t2d)}): mean = {t2d.mean():+.4f}")
    print(f"    ND  (n={len(nd)}):  mean = {nd.mean():+.4f}")
    print(f"    Mann-Whitney U: p = {p:.4e}")
    return score, groups, p, avail


def main():
    # ---- 1. PDX1 probe annotation ----
    print("=" * 74)
    print("1. PDX1 probe annotation on GPL570 (HG-U133 Plus 2.0)")
    print("=" * 74)
    print(f"  PDX1 probes on GPL570: {PDX1_PROBES}")
    for ds in ["GSE76895", "GSE76894"]:
        pres = probe_presence(ds)
        present = [p for p, ok in pres.items() if ok]
        absent = [p for p, ok in pres.items() if not ok]
        print(f"  {ds}: present = {present if present else 'NONE'}; "
              f"absent = {absent}")

    # ---- 2. Uniform 15-gene panel, both cohorts ----
    print()
    print("=" * 74)
    print("2. Uniform 15-gene identity score (16 genes minus PDX1)")
    print("=" * 74)
    score95_15, grp95, p95_15, _ = report("GSE76895", IDG15, "15-gene")
    score94_15, grp94, p94_15, _ = report("GSE76894", IDG15, "15-gene")

    # ---- 3. Discovery cohort: 16-gene vs 15-gene (impact of dropping PDX1) ----
    print()
    print("=" * 74)
    print("3. Discovery cohort (GSE76895): 16-gene vs 15-gene score")
    print("=" * 74)
    score95_16, grp95_16, p95_16, _ = report("GSE76895", IDG16, "16-gene")
    # Pearson correlation between the 15-gene and 16-gene scores.
    r_scores, _ = stats.pearsonr(score95_15, score95_16)
    print(f"  Pearson correlation between 15-gene and 16-gene scores: r = {r_scores:.6f}")

    # ---- Uniform-panel per-sample scores ----
    out_dir = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for ds, score, grp in [("GSE76895", score95_15, grp95),
                           ("GSE76894", score94_15, grp94)]:
        for s, g in zip(score, grp):
            rows.append([ds, g, f"{s:.6f}"])
    out_csv = os.path.join(out_dir, "uniform_15gene_identity_scores.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["dataset", "group", "identity_score_15gene"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
