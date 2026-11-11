# -*- coding: utf-8 -*-
"""
Donor-level single-cell reanalysis for Reviewer 2, Comment 1.

Reviewer 2, Comment 1 (verbatim):
    "Individual cells were treated as independent biological replicates without
     accounting for donor-level variability. Critically, all 495 T2DM beta-cells
     in GSE84133 originated from a single donor, rendering the extremely small
     reported P values potentially misleading. Donor-level pseudobulk or
     mixed-effects analyses and validation in additional diabetic donors are
     essential."

Purpose
-------
This script addresses the comment directly:

1. It documents the donor structure of the Baron (GSE84133) dataset and confirms
   that the 495 T2D beta-cells originate from a single donor (human4). For a
   single-donor group no inferential test is valid, so only a directional log2
   fold change is reported (the cell-level P value reported in the original
   submission is withdrawn).

2. It performs a donor-level pseudobulk re-analysis in the Xin (GSE81608)
   dataset, which provides 6 T2D and 12 non-diabetic donors. Per-donor mean
   expression is computed and compared with a two-tailed Wilcoxon rank-sum test,
   with Benjamini-Hochberg correction across the identity genes. This provides
   the validation in additional diabetic donors requested by the reviewer.

Datasets
--------
1. Baron et al. (2016, Cell Systems; GSE84133)
   - human1 (non-diabetic):  872 beta-cells
   - human2 (non-diabetic):  371 beta-cells
   - human3 (non-diabetic):  787 beta-cells
   - human4 (type 2 diabetic): 495 beta-cells

2. Xin et al. (2016, Cell Metabolism; GSE81608)
   - 12 non-diabetic donors, 6 type 2 diabetic donors (18 donors).

Identity genes
--------------
16 genes (the beta-cell identity panel): PDX1, NKX6-1, HNF1A, HNF4A, MAFB,
NEUROD1, PAX6, SLC2A2, GCK, ABCC8, IAPP, GLP1R, SLC30A8, PCSK1, PCSK2, INS.

Outputs
-------
- baron_donor_structure.txt   : beta-cell count per donor (Baron)
- baron_donor_level_log2fc.csv: donor-level directional log2FC (Baron)
- xin_donor_level_results.csv : donor-level Wilcoxon + BH (Xin)
"""

import csv
import gzip
import os
from collections import defaultdict
from math import erfc, sqrt

import numpy as np

# ---------------------------------------------------------------------------
# Data directory (override with environment variable T2D_SC_DATA if needed).
# By default, walk up from this script to locate the project root containing
# the sc_data directory (consistent with the other response scripts).
# ---------------------------------------------------------------------------
def _find_data_root():
    env = os.environ.get('T2D_SC_DATA')
    if env and os.path.isdir(env):
        return env
    here = os.path.abspath(os.path.dirname(__file__))
    d = here
    for _ in range(12):
        if os.path.isdir(os.path.join(d, 'sc_data')):
            return os.path.join(d, 'sc_data')
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    # Fall back to the old behaviour (script-adjacent sc_data) if not found.
    return os.path.join(here, 'sc_data')

D = _find_data_root()

TARGET = ['PDX1', 'NKX6-1', 'HNF1A', 'HNF4A', 'MAFB', 'NEUROD1', 'PAX6',
          'SLC2A2', 'GCK', 'ABCC8', 'IAPP', 'GLP1R', 'SLC30A8', 'PCSK1', 'PCSK2',
          'INS']


def mannwhitney_u(x, y):
    """Two-sided Mann-Whitney U test using the exact permutation distribution
    when scipy is available, falling back to the asymptotic normal
    approximation with tie correction otherwise. For these sample sizes
    (12 vs 6 donors) the exact p-value is used."""
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
    mu = n1 * n2 / 2.0
    tie_groups = {}
    for v in sorted_vals:
        tie_groups[v] = tie_groups.get(v, 0) + 1
    tie_corr = sum(t ** 3 - t for t in tie_groups.values())
    sigma = sqrt((n1 * n2 / 12.0) * ((n + 1) - tie_corr / (n * (n - 1))))
    z = (U - mu) / sigma
    p = erfc(abs(z) / sqrt(2.0))
    return p


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction (cumulative minimum from the largest
    p-value)."""
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


# ---------------------------------------------------------------------------
# Part 1: Baron (GSE84133) — donor structure and directional log2FC
# ---------------------------------------------------------------------------
def baron_analysis():
    files = {
        'human1': ('ND',  'GSM2230757_human1_umifm_counts.csv.gz'),
        'human2': ('ND',  'GSM2230758_human2_umifm_counts.csv.gz'),
        'human3': ('ND',  'GSM2230759_human3_umifm_counts.csv.gz'),
        'human4': ('T2D', 'GSM2230760_human4_umifm_counts.csv.gz'),
    }

    genes = [g for g in TARGET]
    donor_sum = {g: {} for g in genes}
    donor_ncell = {}

    structure_lines = [
        'Baron et al. (2016, Cell Systems; GSE84133) donor structure:',
    ]

    for tag, (disease, fn) in files.items():
        path = os.path.join(D, fn)
        gene_vals = {g: 0.0 for g in genes}
        ncell = 0
        with gzip.open(path, 'rt') as f:
            header = f.readline().rstrip('\n').split(',')
            gene_idx = {g: header.index(g) for g in genes if g in header}
            cluster_idx = header.index('cluster') if 'cluster' in header else 2
            for line in f:
                parts = line.rstrip('\n').split(',')
                if parts[cluster_idx] == 'beta':
                    ncell += 1
                    for g, gi in gene_idx.items():
                        gene_vals[g] += float(parts[gi])
        donor_ncell[tag] = ncell
        for g in genes:
            donor_sum[g][tag] = gene_vals[g]
        structure_lines.append(f'  {tag} ({disease}): {ncell} beta-cells')

    structure_lines.append('')
    structure_lines.append(
        'The T2D group contains a single donor (human4, 495 beta-cells). '
        'Because n = 1 for the T2D group, no inferential test is valid; only a '
        'directional log2 fold change is reported, and the cell-level P value '
        'from the original submission is withdrawn.'
    )

    out_struct = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'baron_donor_structure.txt')
    with open(out_struct, 'w') as f:
        f.write('\n'.join(structure_lines) + '\n')
    print('\n'.join(structure_lines))
    print()

    out_fc = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'baron_donor_level_log2fc.csv')
    print('=== Baron donor-level log2FC (T2D n=1, directional only) ===')
    with open(out_fc, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'ND_mean_counts', 'T2D_counts', 'log2FC'])
        for g in genes:
            nd_vals = [donor_sum[g][t] for t in ['human1', 'human2', 'human3']]
            nd_mean = np.mean(nd_vals)
            t2d = donor_sum[g]['human4']
            lfc = np.log2((t2d + 0.5) / (nd_mean + 0.5))
            w.writerow([g, f'{nd_mean:.1f}', f'{t2d:.1f}', f'{lfc:.3f}'])
            print(f'{g:10s} {nd_mean:10.1f} {t2d:10.1f} {lfc:+8.3f}')
    print()


# ---------------------------------------------------------------------------
# Part 2: Xin (GSE81608) — donor-level Wilcoxon + BH (multi-donor validation)
# ---------------------------------------------------------------------------
def xin_analysis():
    annot_path = os.path.join(D, 'xin_beta_annot.csv')
    expr_path = os.path.join(D, 'xin_beta_expr.csv')

    # Read annotation: cell -> (donor, condition)
    annot = {}
    with open(annot_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            annot[f'S{row["sample_idx"]}'] = (row['donor'], row['condition'])

    # Read expression: gene -> cell -> RPKM
    expr = {}
    with open(expr_path, encoding='utf-8-sig') as f:
        r = csv.reader(f)
        header = next(r)
        cells = header[1:]
        for row in r:
            gene = row[0]
            expr[gene] = dict(zip(cells, [float(v) for v in row[1:]]))

    donors_nd = sorted({d for d, c in set(annot.values()) if c != 'T2D'})
    donors_t2d = sorted({d for d, c in set(annot.values()) if c == 'T2D'})
    print(f'=== Xin et al. (2016, Cell Metabolism; GSE81608) ===')
    print(f'ND donors ({len(donors_nd)}): {donors_nd}')
    print(f'T2D donors ({len(donors_t2d)}): {donors_t2d}')
    print(f'Total donors: {len(donors_nd) + len(donors_t2d)}')
    print()

    # Donor-level aggregation: mean RPKM per gene per donor
    donor_gene_vals = defaultdict(lambda: defaultdict(list))
    for cell, (donor, cond) in annot.items():
        for gene in TARGET:
            if cell in expr.get(gene, {}):
                donor_gene_vals[donor][gene].append(expr[gene][cell])

    donor_condition = {}
    for cell, (donor, cond) in annot.items():
        donor_condition[donor] = cond

    donor_level = defaultdict(dict)
    for donor, genes in donor_gene_vals.items():
        for gene in TARGET:
            vals = genes.get(gene, [])
            donor_level[donor][gene] = np.mean(vals) if vals else np.nan

    print(f'{"gene":10s} {"ND_mean":>10s} {"T2D_mean":>10s} {"log2FC":>8s} {"p_donor":>10s}')
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

    ps = [r[4] for r in results]
    if ps:
        fdr = bh_fdr(np.array(ps))
        print()
        print('=== Xin donor-level Wilcoxon + BH correction ===')
        for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
            flag = '*' if f < 0.05 else ''
            print(f'{g:10s} log2FC={lfc:+7.3f} p={p:.3e} FDR={f:.3e} {flag}')

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'xin_donor_level_results.csv')
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'ND_mean', 'T2D_mean', 'log2FC', 'p_donor', 'FDR',
                    'ND_n_donor', 'T2D_n_donor'])
        for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
            w.writerow([g, f'{ndm:.4f}', f'{t2dm:.4f}', f'{lfc:.4f}',
                        f'{p:.3e}', f'{f:.3e}', n1, n2])

    print()
    print('Unit of analysis = donor (18 donors: 12 ND + 6 T2D), correcting for cell-level pseudoreplication.')


def main():
    print('#' * 70)
    print('Donor-level single-cell reanalysis (Reviewer 2, Comment 1)')
    print('#' * 70)
    print()
    baron_analysis()
    xin_analysis()


if __name__ == '__main__':
    main()
