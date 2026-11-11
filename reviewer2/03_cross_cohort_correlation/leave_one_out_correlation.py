# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 3:
Cross-cohort correlation, leave-one-gene-out, and robust (rank-based) correlation
analysis for the 15 identity genes common to both bulk cohorts.

This script:
  1. Documents and corrects a methodological error in the original log2FC
     calculation (the expression matrices are already in log2 space after RMA
     normalization, so the correct fold change is the difference of group means
     in log2 space, m_T2D - m_ND, NOT log2(ratio) of log2 values).
  2. Recomputes the corrected per-gene log2FC for GSE76895 and GSE76894.
  3. Computes Pearson and Spearman cross-cohort correlation across all 15 genes.
  4. Performs leave-one-gene-out Pearson correlation.
  5. Performs a direction-consistency (sign concordance) analysis.
"""

import os
import csv
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Locate the data root (the directory containing r3_inputs/).
# The path is resolved from the T2D_DATA environment variable or by walking
# upward from this file.
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(12):
    if os.path.isdir(os.path.join(d, "r3_inputs")):
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
            "Could not locate the data root containing r3_inputs. "
            "Set the T2D_DATA environment variable to the project directory."
        )

R3IN = os.path.join(BASE, "r3_inputs")

# The 16-gene identity panel minus PDX1 (no validation-cohort probe), giving
# the 15 genes common to both platforms.
IDG = [
    "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
    "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8",
    "PCSK1", "PCSK2", "INS",
]


def load_expr(ds):
    """Load the expression matrix (genes x samples) and the sample groups."""
    expr_path = os.path.join(R3IN, f"{ds}_expr.csv")
    grp_path = os.path.join(R3IN, f"{ds}_groups.csv")

    genes = []
    with open(expr_path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])

    gmap = {}
    with open(grp_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gmap[row["sample"]] = row["group"]

    return genes, header[1:], mat, gmap


def corrected_log2fc(ds):
    """Corrected log2 fold change: difference of group means in log2 space.

    The input matrices are RMA-normalized and already in log2 space. Therefore
    the fold change is m_T2D - m_ND (subtraction), not log2(m_T2D / m_ND).
    """
    genes, samples, mat, gmap = load_expr(ds)
    t2d_idx = [j for j, s in enumerate(samples) if gmap.get(s) == "T2D"]
    nd_idx = [j for j, s in enumerate(samples) if gmap.get(s) == "ND"]

    fc = {}
    for g in IDG:
        if g not in genes:
            continue
        gi = genes.index(g)
        m_t = mat[gi, t2d_idx].mean()
        m_n = mat[gi, nd_idx].mean()
        fc[g] = m_t - m_n
    return fc


def main():
    fc95 = corrected_log2fc("GSE76895")
    fc94 = corrected_log2fc("GSE76894")
    common = [g for g in IDG if g in fc95 and g in fc94]

    print("=" * 72)
    print("Corrected per-gene log2FC (m_T2D - m_ND, log2 space)")
    print("=" * 72)
    print(f"{'gene':10s} {'GSE76895':>12s} {'GSE76894':>12s}")
    for g in common:
        print(f"{g:10s} {fc95[g]:+12.4f} {fc94[g]:+12.4f}")

    x = np.array([fc95[g] for g in common])
    y = np.array([fc94[g] for g in common])

    # ---- Pearson / Spearman across all 15 genes ----
    r_pearson, p_pearson = stats.pearsonr(x, y)
    r_spearman, p_spearman = stats.spearmanr(x, y)

    print()
    print("=" * 72)
    print("Cross-cohort correlation (corrected log2FC, 15 genes)")
    print("=" * 72)
    print(f"  Pearson  r = {r_pearson:.4f}  (p = {p_pearson:.3e})")
    print(f"  Spearman rho = {r_spearman:.4f}  (p = {p_spearman:.3e})")

    # ---- Leave-one-gene-out ----
    print()
    print("=" * 72)
    print("Leave-one-gene-out Pearson correlation (corrected log2FC)")
    print("=" * 72)
    loo = {}
    for g in common:
        xl = np.array([fc95[gg] for gg in common if gg != g])
        yl = np.array([fc94[gg] for gg in common if gg != g])
        rl, pl = stats.pearsonr(xl, yl)
        loo[g] = (rl, pl)
        print(f"  excluding {g:10s}: r = {rl:+.4f}  (p = {pl:.3e})")

    # ---- Leave-one-gene-out Spearman (rank-based, robust) ----
    print()
    print("=" * 72)
    print("Leave-one-gene-out Spearman correlation (corrected log2FC)")
    print("=" * 72)
    for g in common:
        xl = np.array([fc95[gg] for gg in common if gg != g])
        yl = np.array([fc94[gg] for gg in common if gg != g])
        rl, pl = stats.spearmanr(xl, yl)
        print(f"  excluding {g:10s}: rho = {rl:+.4f}  (p = {pl:.3e})")

    # ---- Direction-consistency (sign concordance) of the 14 non-SLC2A2 genes ----
    print()
    print("=" * 72)
    print("Direction-consistency of the 14 non-SLC2A2 genes")
    print("=" * 72)
    non_slc = [g for g in common if g != "SLC2A2"]
    concordant = []
    discordant = []
    for g in non_slc:
        s95 = np.sign(fc95[g])
        s94 = np.sign(fc94[g])
        if s95 == s94 and s95 != 0:
            concordant.append(g)
        else:
            discordant.append(g)
    print(f"  Concordant (same direction): {len(concordant)}/{len(non_slc)}")
    print(f"  Discordant (opposite direction): {len(discordant)}/{len(non_slc)}")
    print("  Discordant genes:")
    for g in discordant:
        print(f"    {g:10s}: GSE76895 = {fc95[g]:+.4f}, GSE76894 = {fc94[g]:+.4f}")

    # ---- Write the corrected per-gene log2FC table (Supplementary Table S3) ----
    out_csv = os.path.join(os.path.dirname(__file__),
                           "corrected_log2fc_table.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "log2FC_GSE76895", "log2FC_GSE76894"])
        for g in common:
            w.writerow([g, f"{fc95[g]:.4f}", f"{fc94[g]:.4f}"])


if __name__ == "__main__":
    main()
