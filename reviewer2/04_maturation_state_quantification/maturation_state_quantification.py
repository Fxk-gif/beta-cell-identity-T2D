# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 4:
Direct quantification of beta-cell maturation state and donor-adjusted
proportions, replacing the previous inference from PDX1/NKX6-1 expression.

The reviewer correctly noted that the absence of significant intracellular
changes in PDX1 and NKX6-1 cannot establish depletion of mature beta-cell
subpopulations. We therefore (1) directly quantify a per-cell maturation score
and (2) compare it at the donor level (treating the donor, not the cell, as the
unit of analysis).

Data sources:
  - Xin et al. (GSE81608): purified beta cells, 18 donors (12 ND, 6 T2D),
    503 beta cells. Annotated in sc_data/xin_beta_annot.csv with columns
    sample_idx / condition / donor.
  - Baron et al. (GSE84133): 4 donors (3 ND, 1 T2D). Only one T2D donor, so a
    donor-level statistical comparison is not possible; reported descriptively.

Maturation score definition (per cell):
    maturation score = mean(z-score of mature markers)
                     - mean(z-score of immature/progenitor markers)

  Mature markers:      MAFA, UCN3, INS, SLC2A2, GCK, SLC30A8, PCSK1
  Immature markers:    MAFB, NEUROG3, CHL1, MYCN, OPRD1

z-scores are computed per gene across all beta cells in the dataset.
"""

import os
import csv
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Locate the data root (the directory containing sc_data/) by walking upward
# from this file.
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(12):
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

# Maturation marker sets. Only genes present in the expression matrix are used.
MATURE = ["MAFA", "UCN3", "INS", "SLC2A2", "GCK", "SLC30A8", "PCSK1"]
IMMATURE = ["MAFB", "NEUROG3", "CHL1", "MYCN", "OPRD1"]


def load_xin():
    """Load Xin beta-cell expression and donor annotation.

    Returns (genes, cell_ids, mat, donor_of_cell, cond_of_cell).
    The expression matrix is genes x cells; cells are S0..S(n-1).
    """
    expr_path = os.path.join(SCD, "xin_beta_expr.csv")
    annot_path = os.path.join(SCD, "xin_beta_annot.csv")

    with open(expr_path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    cell_ids = header[1:]  # e.g. 'S0', 'S1', ...

    # Annotation maps sample_idx (int, 0-based) -> (donor, condition).
    donor_of_cell = []
    cond_of_cell = []
    idx_to_donor = {}
    idx_to_cond = {}
    with open(annot_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            i = int(row["sample_idx"])
            idx_to_donor[i] = row["donor"]
            idx_to_cond[i] = row["condition"]
    # cell_ids are S<idx>; map them back to the integer index.
    for cid in cell_ids:
        i = int(cid.lstrip("S"))
        donor_of_cell.append(idx_to_donor.get(i, "NA"))
        cond_of_cell.append(idx_to_cond.get(i, "NA"))

    return genes, cell_ids, mat, donor_of_cell, cond_of_cell


def zscore_rows(mat):
    """z-score each row (gene) across all columns (cells), ddof=1.

    Genes with zero variance across cells (e.g. a marker that is not expressed
    in any mature beta cell, such as NEUROG3) receive a z-score of 0, so that
    they do not contribute NaN to the maturation score.
    """
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True, ddof=1)
    z = np.zeros_like(mat, dtype=float)
    ok = (sd[:, 0] > 0)
    z[ok] = (mat[ok] - mu[ok]) / sd[ok]
    return z


def maturation_score(mat, genes):
    """Per-cell maturation score from the expression matrix (genes x cells)."""
    z = zscore_rows(mat)
    gi = {g: i for i, g in enumerate(genes)}

    mat_genes = [g for g in MATURE if g in gi]
    imm_genes = [g for g in IMMATURE if g in gi]

    mat_idx = [gi[g] for g in mat_genes]
    imm_idx = [gi[g] for g in imm_genes]

    mat_mean = z[mat_idx].mean(axis=0)
    imm_mean = z[imm_idx].mean(axis=0)
    score = mat_mean - imm_mean
    return score, mat_genes, imm_genes


def donor_level(score, donor_of_cell, cond_of_cell):
    """Aggregate per-cell score to per-donor mean score."""
    donors = sorted(set(donor_of_cell))
    rows = []
    for dn in donors:
        idx = [i for i, dd in enumerate(donor_of_cell) if dd == dn]
        s = score[idx]
        cond = cond_of_cell[idx[0]]
        rows.append({
            "donor": dn,
            "condition": cond,
            "n_cells": len(idx),
            "mean_score": float(s.mean()),
        })
    return rows


def main():
    genes, cell_ids, mat, donor_of_cell, cond_of_cell = load_xin()

    score, mat_genes, imm_genes = maturation_score(mat, genes)

    print("=" * 74)
    print("Direct quantification of beta-cell maturation state (Xin, GSE81608)")
    print("=" * 74)
    print(f"  Mature markers used:   {mat_genes}")
    print(f"  Immature markers used: {imm_genes}")
    print(f"  Total beta cells:      {len(cell_ids)}")

    # ---- Per-cell maturation score: T2D vs ND (cell-level, descriptive) ----
    cond_arr = np.array(cond_of_cell)
    t2d_mask = cond_arr == "T2D"
    nd_mask = cond_arr == "non-diabetic"
    t2d_score = score[t2d_mask]
    nd_score = score[nd_mask]

    # Mann-Whitney U at cell level (for completeness; donor-level is primary)
    _, p_cell = stats.mannwhitneyu(t2d_score, nd_score, alternative="two-sided")

    print()
    print("Cell-level maturation score (descriptive, not the primary inference):")
    print(f"  T2D beta cells (n={t2d_mask.sum()}): mean = {t2d_score.mean():+.4f}")
    print(f"  ND  beta cells (n={nd_mask.sum()}):  mean = {nd_score.mean():+.4f}")
    print(f"  Mann-Whitney U (cell-level): p = {p_cell:.4f}")

    # ---- Donor-level maturation score (PRIMARY inference) ----
    dl = donor_level(score, donor_of_cell, cond_of_cell)

    print()
    print("=" * 74)
    print("Donor-level mean maturation score (primary analysis)")
    print("=" * 74)
    t2d_donors = [r for r in dl if r["condition"] == "T2D"]
    nd_donors = [r for r in dl if r["condition"] == "non-diabetic"]

    for r in dl:
        print(f"  {r['donor']:12s} ({r['condition']:14s}): "
              f"n={r['n_cells']:3d}, mean_score={r['mean_score']:+.4f}")

    t2d_vals = np.array([r["mean_score"] for r in t2d_donors])
    nd_vals = np.array([r["mean_score"] for r in nd_donors])

    _, p_d = stats.mannwhitneyu(t2d_vals, nd_vals, alternative="two-sided")

    print()
    print(f"  T2D donors (n={len(t2d_vals)}): mean = {t2d_vals.mean():+.4f}, "
          f"sd = {t2d_vals.std(ddof=1):.4f}")
    print(f"  ND  donors (n={len(nd_vals)}):  mean = {nd_vals.mean():+.4f}, "
          f"sd = {nd_vals.std(ddof=1):.4f}")
    print(f"  Mann-Whitney U (donor-level): p = {p_d:.4f}")

    # ---- Donor-level mature beta-cell proportion ----
    # Classify a cell as "mature" if its score is above the median score of
    # ND cells (an ND-derived threshold); then compute per-donor proportion.
    nd_threshold = float(np.median(nd_score))
    is_mature = score > nd_threshold

    print()
    print("=" * 74)
    print("Donor-level mature beta-cell proportion (threshold = ND median)")
    print("=" * 74)
    print(f"  ND-derived maturity threshold (median ND score) = {nd_threshold:+.4f}")
    prop_rows = []
    for dn in sorted(set(donor_of_cell)):
        idx = [i for i, dd in enumerate(donor_of_cell) if dd == dn]
        pct = 100.0 * is_mature[idx].mean()
        cond = cond_of_cell[idx[0]]
        prop_rows.append({"donor": dn, "condition": cond, "mature_pct": pct})
        print(f"  {dn:12s} ({cond:14s}): mature beta fraction = {pct:5.2f}%")

    t2d_prop = np.array([r["mature_pct"] for r in prop_rows if r["condition"] == "T2D"])
    nd_prop = np.array([r["mature_pct"] for r in prop_rows if r["condition"] == "non-diabetic"])
    _, p_prop = stats.mannwhitneyu(t2d_prop, nd_prop, alternative="two-sided")

    print()
    print(f"  T2D donor mature proportion: mean = {t2d_prop.mean():.2f}%")
    print(f"  ND  donor mature proportion: mean = {nd_prop.mean():.2f}%")
    print(f"  Mann-Whitney U (donor-level): p = {p_prop:.4f}")

    # ---- Individual maturation markers (donor-level) ----
    print()
    print("=" * 74)
    print("Individual maturation markers (cell-level, T2D vs ND)")
    print("=" * 74)
    for g in ["MAFA", "MAFB", "UCN3", "NEUROG3", "CHL1"]:
        if g not in genes:
            print(f"  {g}: not present in matrix")
            continue
        gi = genes.index(g)
        v = mat[gi]
        vt = v[t2d_mask]
        vn = v[nd_mask]
        # log2 fold change (expression matrix appears to be in linear/RPKM-like
        # space for Xin; compute log2 ratio of means with a small offset)
        fc = np.log2((vt.mean() + 0.01) / (vn.mean() + 0.01))
        _, p_g = stats.mannwhitneyu(vt, vn, alternative="two-sided")
        print(f"  {g:10s}: T2D_mean={vt.mean():.4f}, ND_mean={vn.mean():.4f}, "
              f"log2FC={fc:+.3f}, p={p_g:.4f}")

    # ---- Per-donor result tables ----
    out_csv = os.path.join(os.path.dirname(__file__),
                           "maturation_score_donor_level.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["donor", "condition", "n_cells", "mean_maturation_score"])
        for r in dl:
            w.writerow([r["donor"], r["condition"], r["n_cells"],
                        f"{r['mean_score']:.6f}"])

    out_prop = os.path.join(os.path.dirname(__file__),
                            "mature_beta_proportion_donor_level.csv")
    with open(out_prop, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["donor", "condition", "mature_beta_pct"])
        for r in prop_rows:
            w.writerow([r["donor"], r["condition"], f"{r['mature_pct']:.4f}"])


if __name__ == "__main__":
    main()
