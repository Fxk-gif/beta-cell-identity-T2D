# -*- coding: utf-8 -*-
"""
Bulk progenitor/disallowed score and the trace-level caveat.

The eight genes in this panel are expressed at trace levels in bulk islets, so
a highly significant p-value on near-floor microarray signal is weak evidence.
This script recomputes the bulk progenitor/disallowed score at full precision
and reports the single-gene contributions that underlie it, so that the trace
level of each gene can be inspected directly.

Score definition:
    progenitor/disallowed score = mean(z-score of the 8 genes), where the
    z-score of each gene is computed across samples within a cohort.

    Gene panel (8 canonical fetal/progenitor and disallowed genes):
        FOXE1, MYCN, CHL1, OPRD1, GRIN2D, ALDH1A3, HES1, SOX9

    Only genes with a probe on the GPL570 platform are retained. FOXE1 has no
    probe in the validation cohort (GSE76894), so the validation score is
    computed from the remaining 7 genes; this is stated in the output.

Statistical test:
    Wilcoxon rank-sum (Mann-Whitney U) test, exact method (matching
    wilcox.test(..., exact = TRUE) in R).

Interpretation (see response letter, Comment 8):
    The score is lower in T2D than in non-diabetic islets in both cohorts, a
    direction opposite to an active dedifferentiation program. Because the
    constituent genes are expressed at trace levels, this reduction is read as
    reflecting the same depletion of mature beta cells rather than as strong
    evidence against dedifferentiation.
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

R3 = os.path.join(BASE, "r3_inputs")

# The 8-gene progenitor/disallowed panel.
PROGENITOR_GENES = ["FOXE1", "MYCN", "CHL1", "OPRD1", "GRIN2D",
                    "ALDH1A3", "HES1", "SOX9"]

COHORTS = [
    ("GSE76895", "GSE76895_expr.csv", "GSE76895_groups.csv"),
    ("GSE76894", "GSE76894_expr.csv", "GSE76894_groups.csv"),
]


def zscore(x):
    """Per-gene z-score across samples (n-1 standard deviation)."""
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)


def load_cohort(expr_file, grp_file):
    """Load a bulk cohort, returning (genes, sample_ids, matrix, groups).

    The expression matrix is genes x samples; groups is a dict sample -> group.
    """
    genes = []
    samples = []
    rows = []
    with open(os.path.join(R3, expr_file), newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        for line in reader:
            genes.append(line[0])
            rows.append([float(v) for v in line[1:]])
    mat = np.asarray(rows)

    groups = {}
    with open(os.path.join(R3, grp_file), newline="") as f:
        reader = csv.reader(f)
        next(reader)  # header
        for line in reader:
            if line:
                groups[line[0]] = line[1]
    return genes, samples, mat, groups


def main():
    results = []
    cohort_p = {}
    out_lines = []
    out_lines.append("Bulk progenitor/disallowed score (8-gene panel)")
    out_lines.append("Panel: " + ", ".join(PROGENITOR_GENES))
    out_lines.append("Test: Wilcoxon rank-sum (exact), T2D vs ND")
    out_lines.append("")

    for tag, expr_file, grp_file in COHORTS:
        genes, samples, mat, groups = load_cohort(expr_file, grp_file)
        gene_idx = {g: i for i, g in enumerate(genes)}

        avail = [g for g in PROGENITOR_GENES if g in gene_idx]
        missing = [g for g in PROGENITOR_GENES if g not in gene_idx]

        keep = [i for i, s in enumerate(samples) if groups.get(s) in ("T2D", "ND")]
        samples_k = [samples[i] for i in keep]
        mat_k = mat[:, keep]
        grp_arr = np.array([1 if groups[s] == "T2D" else 0 for s in samples_k])

        # z-score each available gene across samples, then average across genes.
        sub = np.vstack([zscore(mat_k[gene_idx[g]]) for g in avail])
        score = sub.mean(axis=0)

        t2d = grp_arr == 1
        nd = grp_arr == 0
        p_exact = stats.mannwhitneyu(
            score[t2d], score[nd], alternative="two-sided", method="exact"
        ).pvalue
        p_asymp = stats.mannwhitneyu(
            score[t2d], score[nd], alternative="two-sided", method="asymptotic"
        ).pvalue
        cohort_p[tag] = p_exact

        out_lines.append("=" * 60)
        out_lines.append(f"Cohort: {tag}  (T2D n={int(t2d.sum())}, ND n={int(nd.sum())})")
        out_lines.append(f"Genes used ({len(avail)}): " + ", ".join(avail))
        if missing:
            out_lines.append("Genes absent from platform: " + ", ".join(missing))
        out_lines.append(f"Score mean  ND = {score[nd].mean():+.4f}, "
                         f"T2D = {score[t2d].mean():+.4f}")
        out_lines.append(f"Wilcoxon exact p = {p_exact:.6e}")
        out_lines.append(f"Wilcoxon asymptotic p = {p_asymp:.6e} (reference)")
        out_lines.append("")
        out_lines.append("Single-gene contribution (group-mean expression, exact p):")
        out_lines.append(f"{'gene':10s} {'ND mean':>10s} {'T2D mean':>10s} {'exact p':>12s}")
        for g in avail:
            vals = mat_k[gene_idx[g]]
            gp = stats.mannwhitneyu(
                vals[t2d], vals[nd], alternative="two-sided", method="exact"
            ).pvalue
            out_lines.append(f"{g:10s} {vals[nd].mean():>10.3f} "
                             f"{vals[t2d].mean():>10.3f} {gp:>12.3e}")
            results.append([tag, g, f"{vals[nd].mean():.3f}",
                            f"{vals[t2d].mean():.3f}", f"{gp:.6e}"])
        # score-level summary row
        results.append([tag, "__SCORE__", f"{score[nd].mean():.4f}",
                        f"{score[t2d].mean():.4f}", f"{p_exact:.6e}"])
        out_lines.append("")

    out_lines.append("=" * 60)
    out_lines.append("Summary (full-precision p-values, exact method):")
    for tag, _expr, _grp in COHORTS:
        out_lines.append(f"  {tag}: p = {cohort_p[tag]:.3e}")
    out_lines.append("")
    out_lines.append("Note: FOXE1 is absent from the GPL570 probes of GSE76894, so")
    out_lines.append("the validation score is computed from the remaining 7 genes.")

    text = "\n".join(out_lines)
    out_path = os.path.join(os.path.dirname(__file__), "progenitor_score_output.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)

    csv_path = os.path.join(os.path.dirname(__file__), "progenitor_score_results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cohort", "gene", "ND_mean", "T2D_mean", "exact_p"])
        w.writerows(results)


if __name__ == "__main__":
    main()
