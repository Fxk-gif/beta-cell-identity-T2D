# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 7:
Verification of the full identity score in the Xin single-cell dataset, and
clarification that GSEA/PPI establish association rather than a causal cascade.

The reviewer stated that the full identity score was nonsignificant in the Xin
single-cell dataset. The full 16-gene identity score (mean of per-gene
z-scores) is recomputed here and tested for T2D vs ND at both the cell level
and the donor level, reporting the P value and effect size at each level.

Identity score definition (consistent with the bulk cohorts):
    identity score = mean( z-score of each of the 16 identity genes )
    z-scores are computed per gene across all beta cells in the dataset.

The 16 identity genes are:
    PDX1, NKX6-1, HNF1A, HNF4A, MAFB, NEUROD1, PAX6,
    SLC2A2, GCK, ABCC8, IAPP, GLP1R, SLC30A8, PCSK1, PCSK2, INS
"""

import os
import csv
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Locate the data root (directory containing sc_data).
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(10):
    if os.path.isdir(os.path.join(d, "sc_data")):
        BASE = d
        break
    parent = os.path.dirname(d)
    if parent == d:
        break
    d = parent

if BASE is None:
    BASE = os.environ.get("T2D_DATA", "")
    if not BASE or not os.path.isdir(os.path.join(BASE, "sc_data")):
        raise FileNotFoundError(
            "Could not locate the data root containing sc_data. "
            "Set the T2D_DATA environment variable to the project directory."
        )

SCD = os.path.join(BASE, "sc_data")

# 16-gene identity panel (the "full identity score").
IDG16 = [
    "PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
    "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8", "PCSK1",
    "PCSK2", "INS",
]


def load_xin():
    expr_path = os.path.join(SCD, "xin_beta_expr.csv")
    annot_path = os.path.join(SCD, "xin_beta_annot.csv")

    with open(expr_path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    cell_ids = header[1:]

    donor_of_cell = []
    cond_of_cell = []
    idx_to_donor = {}
    idx_to_cond = {}
    with open(annot_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            i = int(row["sample_idx"])
            idx_to_donor[i] = row["donor"]
            idx_to_cond[i] = row["condition"]
    for cid in cell_ids:
        i = int(cid.lstrip("S"))
        donor_of_cell.append(idx_to_donor.get(i, "NA"))
        cond_of_cell.append(idx_to_cond.get(i, "NA"))

    return genes, cell_ids, mat, donor_of_cell, cond_of_cell


def zscore_rows(mat):
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True, ddof=1)
    z = np.zeros_like(mat, dtype=float)
    ok = (sd[:, 0] > 0)
    z[ok] = (mat[ok] - mu[ok]) / sd[ok]
    return z


def identity_score(mat, genes):
    z = zscore_rows(mat)
    gi = {g: i for i, g in enumerate(genes)}
    used = [g for g in IDG16 if g in gi]
    idx = [gi[g] for g in used]
    score = z[idx].mean(axis=0)
    return score, used


def main():
    genes, cell_ids, mat, donor_of_cell, cond_of_cell = load_xin()

    score, used = identity_score(mat, genes)
    missing = [g for g in IDG16 if g not in genes]

    print("=" * 76)
    print("Full identity score in the Xin single-cell dataset (GSE81608)")
    print("=" * 76)
    print(f"  Identity genes used:  {len(used)} / {len(IDG16)}")
    print(f"  Genes used:           {used}")
    print(f"  Genes missing:        {missing if missing else 'none'}")
    print(f"  Total beta cells:     {len(cell_ids)}")

    cond_arr = np.array(cond_of_cell)
    t2d_mask = cond_arr == "T2D"
    nd_mask = cond_arr == "non-diabetic"
    t2d_score = score[t2d_mask]
    nd_score = score[nd_mask]

    # ---- Cell-level (descriptive; donor-level is primary) ----
    u_cell, p_cell = stats.mannwhitneyu(t2d_score, nd_score, alternative="two-sided")
    t_cell, p_t_cell = stats.ttest_ind(t2d_score, nd_score, equal_var=False)

    print()
    print("Cell-level identity score (descriptive):")
    print(f"  T2D beta cells (n={t2d_mask.sum()}): mean = {t2d_score.mean():+.4f}, "
          f"sd = {t2d_score.std(ddof=1):.4f}")
    print(f"  ND  beta cells (n={nd_mask.sum()}):  mean = {nd_score.mean():+.4f}, "
          f"sd = {nd_score.std(ddof=1):.4f}")
    print(f"  Mann-Whitney U (cell-level): p = {p_cell:.4f}")
    print(f"  Welch t-test    (cell-level): p = {p_t_cell:.4f}")

    # ---- Donor-level (primary inference; treats donor as unit) ----
    donors = sorted(set(donor_of_cell))
    dl = []
    for dn in donors:
        idx = [i for i, dd in enumerate(donor_of_cell) if dd == dn]
        dl.append({
            "donor": dn,
            "condition": cond_of_cell[idx[0]],
            "n_cells": len(idx),
            "mean_score": float(score[idx].mean()),
        })

    t2d_donors = [r for r in dl if r["condition"] == "T2D"]
    nd_donors = [r for r in dl if r["condition"] == "non-diabetic"]

    print()
    print("=" * 76)
    print("Donor-level mean identity score (primary analysis)")
    print("=" * 76)
    for r in dl:
        print(f"  {r['donor']:12s} ({r['condition']:14s}): "
              f"n={r['n_cells']:3d}, mean_score={r['mean_score']:+.4f}")

    t2d_vals = np.array([r["mean_score"] for r in t2d_donors])
    nd_vals = np.array([r["mean_score"] for r in nd_donors])

    u_d, p_d = stats.mannwhitneyu(t2d_vals, nd_vals, alternative="two-sided")
    t_d, p_t_d = stats.ttest_ind(t2d_vals, nd_vals, equal_var=False)

    print()
    print(f"  T2D donors (n={len(t2d_vals)}): mean = {t2d_vals.mean():+.4f}, "
          f"sd = {t2d_vals.std(ddof=1):.4f}")
    print(f"  ND  donors (n={len(nd_vals)}):  mean = {nd_vals.mean():+.4f}, "
          f"sd = {nd_vals.std(ddof=1):.4f}")
    print(f"  Mann-Whitney U (donor-level): p = {p_d:.4f}")
    print(f"  Welch t-test    (donor-level): p = {p_t_d:.4f}")

    # Effect size (Cohen's d) at donor level.
    pooled_sd = np.sqrt(((len(t2d_vals) - 1) * t2d_vals.std(ddof=1) ** 2 +
                         (len(nd_vals) - 1) * nd_vals.std(ddof=1) ** 2) /
                        (len(t2d_vals) + len(nd_vals) - 2))
    cohen_d = (t2d_vals.mean() - nd_vals.mean()) / pooled_sd if pooled_sd > 0 else float("nan")
    print(f"  Cohen's d (donor-level): {cohen_d:+.4f}")

    # ---- Individual genes at donor level (which genes drive any signal) ----
    # Per-donor means of the raw RPKM expression (pseudobulk), followed by a
    # donor-level Mann-Whitney U test. scipy method="auto" uses the exact
    # distribution for small samples.
    print()
    print("=" * 76)
    print("Donor-level individual genes (raw-RPKM pseudobulk, Mann-Whitney U)")
    print("=" * 76)
    gi = {g: i for i, g in enumerate(genes)}
    gene_results = []
    for g in used:
        idx = gi[g]
        donor_mean_t = []
        donor_mean_n = []
        for r in dl:
            ci = [i for i, dd in enumerate(donor_of_cell) if dd == r["donor"]]
            v = mat[idx, ci].mean()  # per-donor mean of raw RPKM (not z-score)
            if r["condition"] == "T2D":
                donor_mean_t.append(v)
            else:
                donor_mean_n.append(v)
        u_g, p_g = stats.mannwhitneyu(donor_mean_t, donor_mean_n,
                                      alternative="two-sided")
        # Donor-level log2FC (T2D mean − ND mean, in log2 space)
        log2fc = np.log2(np.mean(donor_mean_t) / np.mean(donor_mean_n))
        gene_results.append((g, log2fc, p_g))

    gene_results.sort(key=lambda x: x[2])
    print(f"{'gene':10s} {'log2FC':>10s} {'p(donor-level)':>15s}")
    for g, log2fc, p_g in gene_results:
        print(f"{g:10s} {log2fc:+10.4f} {p_g:15.4f}")

    # ---- Donor-level identity scores (CSV output) ----
    out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "xin_identity_score_donor_level.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["donor", "condition", "n_cells", "mean_identity_score"])
        for r in dl:
            w.writerow([r["donor"], r["condition"], r["n_cells"],
                        f"{r['mean_score']:.6f}"])

    print()


if __name__ == "__main__":
    main()
