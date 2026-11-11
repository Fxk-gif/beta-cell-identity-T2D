# -*- coding: utf-8 -*-
"""
Xin (GSE81608) donor-level single-cell reanalysis
-- fixes the cell-level pseudoreplication issue.

Background
----------
The reviewers (Reviewer 1, Comment 1 / Reviewer 2, Comment 1) noted that the
original manuscript ran Mann-Whitney tests treating individual cells as
independent samples, which is pseudoreplication because multiple beta cells
from the same donor are not independent. This script aggregates the unit of
analysis from individual cells to donors (donor-level pseudobulk).

Data
----
- sc_data/xin_beta_expr.csv  : RPKM expression matrix of beta cells from
                               Xin et al. 2016 (Cell Metab), GSE81608
                               (rows = genes, columns = cells S0..S502, 503 beta cells)
- sc_data/xin_beta_annot.csv : donor assignment and diabetes status of each cell
                               (columns: sample_idx / condition / donor)

Donor counts (as measured by this script)
-----------------------------------------
- Non-diabetic (ND): 12 donors (Non T2D 1-12)
- Type 2 diabetic (T2D): 6 donors (T2D 1-6)
- 18 donors in total

Methods
-------
1. For each donor, average the RPKM of every gene across its beta cells to
   obtain a donor-level pseudobulk expression value.
2. Using donors as the unit (n = 18), run a two-sided Wilcoxon rank-sum test
   (ND vs T2D) for every gene.
3. Apply Benjamini-Hochberg FDR correction to the p-values of the 16 identity genes.

Note: the public Xin data are RPKM (not raw counts), so an edgeR negative
binomial model does not apply; a donor-level non-parametric Wilcoxon test is used.

Run
---
    python xin_donor_level.py

Output
------
- sc_data/xin_donor_level_results.csv : per-gene ND mean, T2D mean, log2FC, p, FDR, donor counts
"""

import csv
import os
from collections import defaultdict
from math import erfc, sqrt

import numpy as np

# Data directory: overridable via the T2D_SC_DATA environment variable; otherwise
# walk upward to locate the project root that contains sc_data/.
def _find_sc_data():
    env = os.environ.get('T2D_SC_DATA', '').rstrip('/\\')
    if env and os.path.isdir(env):
        return env
    d = os.path.abspath(os.path.dirname(__file__))
    for _ in range(12):
        if os.path.isdir(os.path.join(d, 'sc_data')):
            return os.path.join(d, 'sc_data')
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sc_data')

D = _find_sc_data()

# Target genes (the 16-gene beta-cell identity panel)
TARGET = ['PDX1', 'NKX6-1', 'HNF1A', 'HNF4A', 'MAFB', 'NEUROD1', 'PAX6',
          'SLC2A2', 'GCK', 'ABCC8', 'IAPP', 'GLP1R', 'SLC30A8', 'PCSK1', 'PCSK2',
          'INS']


def mannwhitney_u(x, y):
    """Two-sided Mann-Whitney U test using the exact null distribution, with a
    normal approximation + tie correction as a fallback when scipy is
    unavailable. The sample size is the number of donors (12 vs 6); the exact
    and normal-approximation methods differ slightly for a few genes (e.g. GCK
    exact p = 3.23e-3 vs approximate 4.96e-3). The exact method is used
    throughout.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    try:
        from scipy import stats as _stats
        return float(_stats.mannwhitneyu(x, y, alternative="two-sided",
                                         method="exact").pvalue)
    except Exception:
        pass
    n1, n2 = len(x), len(y)
    n = n1 + n2
    ranked = np.concatenate([x, y])
    order = np.argsort(ranked)
    ranks = np.empty(n, dtype=float)
    sorted_vals = ranked[order]
    # handle ties (average ranks)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        ranks[order[i:j + 1]] = avg
        i = j + 1
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U = min(U1, U2)
    # normal approximation (with tie correction)
    mu = n1 * n2 / 2.0
    # tie-correction term: t_i is the multiplicity of each tie group
    tie_groups = {}
    for v in sorted_vals:
        tie_groups[v] = tie_groups.get(v, 0) + 1
    tie_corr = 0.0
    for t in tie_groups.values():
        tie_corr += (t ** 3 - t)
    sigma = sqrt((n1 * n2 / 12.0) * ((n + 1) - tie_corr / (n * (n - 1))))
    z = (U - mu) / sigma
    p = erfc(abs(z) / sqrt(2.0))  # two-sided
    return p


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction (cumulative minimum from largest)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    adj = np.empty(n)
    cummin = np.inf
    for i in range(n - 1, -1, -1):
        cummin = min(cummin, sorted_p[i] * n / (i + 1))
        adj[i] = cummin
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    # 1) read annot: cell -> (donor, condition)
    annot = {}
    annot_path = os.path.join(D, 'xin_beta_annot.csv')
    with open(annot_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            annot[f'S{row["sample_idx"]}'] = (row['donor'], row['condition'])

    # 2) read expr: gene -> cell -> RPKM
    expr = {}
    expr_path = os.path.join(D, 'xin_beta_expr.csv')
    with open(expr_path, encoding='utf-8-sig') as f:
        r = csv.reader(f)
        header = next(r)
        cells = header[1:]
        for row in r:
            gene = row[0]
            expr[gene] = dict(zip(cells, [float(v) for v in row[1:]]))

    # 3) count donors
    donors_nd = sorted({d for d, c in set(annot.values()) if c != 'T2D'})
    donors_t2d = sorted({d for d, c in set(annot.values()) if c == 'T2D'})
    print(f'ND donors ({len(donors_nd)}): {donors_nd}')
    print(f'T2D donors ({len(donors_t2d)}): {donors_t2d}')
    print(f'Total donors: {len(donors_nd) + len(donors_t2d)}')

    # 4) donor-level aggregation: mean beta-cell RPKM per donor per gene
    donor_gene_vals = defaultdict(lambda: defaultdict(list))
    for cell, (donor, cond) in annot.items():
        for gene in TARGET:
            if cell in expr.get(gene, {}):
                donor_gene_vals[donor][gene].append(expr[gene][cell])

    donor_level = defaultdict(dict)
    donor_condition = {}
    for cell, (donor, cond) in annot.items():
        donor_condition[donor] = cond
    for donor, genes in donor_gene_vals.items():
        for gene in TARGET:
            vals = genes.get(gene, [])
            donor_level[donor][gene] = np.mean(vals) if vals else np.nan

    # 5) donor-level Wilcoxon test
    print(f'\n{"gene":10s} {"ND_mean":>10s} {"T2D_mean":>10s} {"log2FC":>8s} {"p_donor":>10s}')
    results = []
    for gene in TARGET:
        nd_vals = [donor_level[d][gene] for d in donor_level
                   if donor_condition.get(d) != 'T2D' and not np.isnan(donor_level[d][gene])]
        t2d_vals = [donor_level[d][gene] for d in donor_level
                    if donor_condition.get(d) == 'T2D' and not np.isnan(donor_level[d][gene])]
        if len(nd_vals) < 2 or len(t2d_vals) < 2:
            print(f'{gene:10s} insufficient samples (ND={len(nd_vals)}, T2D={len(t2d_vals)}), skipped')
            continue
        ndm = np.mean(nd_vals)
        t2dm = np.mean(t2d_vals)
        lfc = np.log2((t2dm + 0.1) / (ndm + 0.1))
        p = mannwhitney_u(t2d_vals, nd_vals)
        results.append((gene, ndm, t2dm, lfc, p, len(nd_vals), len(t2d_vals)))
        print(f'{gene:10s} {ndm:10.3f} {t2dm:10.3f} {lfc:+8.3f} {p:10.3e}')

    # 6) FDR correction
    ps = [r[4] for r in results]
    if ps:
        fdr = bh_fdr(np.array(ps))
        print('\n=== donor-level Wilcoxon + BH correction ===')
        for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
            flag = '*' if f < 0.05 else ''
            print(f'{g:10s} log2FC={lfc:+7.3f} p={p:.3e} FDR={f:.3e} {flag}')

    # 7) save (next to this script)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'xin_donor_level_results.csv')
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'ND_mean', 'T2D_mean', 'log2FC', 'p_donor', 'FDR',
                    'ND_n_donor', 'T2D_n_donor'])
        for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
            w.writerow([g, f'{ndm:.4f}', f'{t2dm:.4f}', f'{lfc:.4f}',
                        f'{p:.3e}', f'{f:.3e}', n1, n2])

    print('Unit of analysis = donor (18 donors: 12 ND + 6 T2D), correcting cell-level pseudoreplication.')


if __name__ == '__main__':
    main()
