# -*- coding: utf-8 -*-
"""
Baron (GSE84133) donor-level single-cell reanalysis
-- fixes the cell-level pseudoreplication issue.

Background
----------
The reviewers (Reviewer 1, Comment 1 / Reviewer 2, Comment 1) noted that the
T2DM group in the Baron dataset consists of a single donor, yet the original
manuscript ran a Mann-Whitney test across cells and reported FDR = 2.4e-34,
which is pseudoreplication. With n = 1 donor, the only defensible statistic is
a fold change.

This script retains donor identity and aggregates beta cells per donor. Because
the T2D group has only 1 donor, it reports only a directional log2FC and runs
no statistical test.

Data
----
- sc_data/GSM2230757_human1_umifm_counts.csv.gz  (human1, non-diabetic)
- sc_data/GSM2230758_human2_umifm_counts.csv.gz  (human2, non-diabetic)
- sc_data/GSM2230759_human3_umifm_counts.csv.gz  (human3, non-diabetic)
- sc_data/GSM2230760_human4_umifm_counts.csv.gz  (human4, type 2 diabetic)

(Baron et al. 2016, Cell Systems; GSE84133)

Donor and beta-cell counts (as measured by this script)
-------------------------------------------------------
- human1 (ND) : 872 beta cells
- human2 (ND) : 371 beta cells
- human3 (ND) : 787 beta cells
- human4 (T2D): 495 beta cells
  -> the T2D group has only 1 donor (human4), so no statistical significance is
     reported; only a directional log2FC is given.

Methods
-------
1. Select beta cells (cluster == 'beta') from each donor's raw counts file.
2. Sum counts per donor per gene to obtain donor-level pseudobulk counts.
3. Compute log2FC = log2( (T2D_counts + 0.5) / (ND_mean_counts + 0.5) ).
4. Do not report p-values / FDR (because T2D n = 1).

Run
---
    python baron_donor_level.py

Output
------
- sc_data/baron_donor_level_counts.csv : total beta-cell counts per gene across the 4 donors
- sc_data/baron_donor_level_log2fc.csv : per-gene ND mean, T2D counts, log2FC
"""

import csv
import gzip
import os

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

# Donor identity and disease status, and source files
DONOR_DISEASE = {'human1': 'ND', 'human2': 'ND', 'human3': 'ND', 'human4': 'T2D'}
FILES = {
    'human1': 'GSM2230757_human1_umifm_counts.csv.gz',
    'human2': 'GSM2230758_human2_umifm_counts.csv.gz',
    'human3': 'GSM2230759_human3_umifm_counts.csv.gz',
    'human4': 'GSM2230760_human4_umifm_counts.csv.gz',
}

# Identity genes (16 core identity genes) plus 8
# progenitor/disallowed-related genes detectable in the Baron single-cell data
# (single-cell trace-level expression probes; note that these 8 are the subset
# measurable in the single-cell data, and differ in coverage from the bulk-level
# 8-gene score FOXE1/MYCN/CHL1/OPRD1/GRIN2D/ALDH1A3/HES1/SOX9; they are used
# here only as a single-cell directional comparison)
TARGET = ['SLC2A2', 'IAPP', 'PDX1', 'NKX6-1', 'HNF1A', 'HNF4A', 'MAFA', 'MAFB', 'NEUROD1',
          'GCK', 'GLP1R', 'INS', 'ABCC8', 'KCNJ11', 'SLC30A8', 'PCSK2', 'CHL1',
          'GAD1', 'CALCA', 'TSHR', 'NHLH2', 'FOXE1', 'OPRD1', 'GRIN2D']


def main():
    donor_beta_counts = {}
    donor_beta_ncell = {}

    for tag, fn in FILES.items():
        path = os.path.join(D, fn)
        gene_lists = {g: [] for g in TARGET}
        ncell = 0
        with gzip.open(path, 'rt') as f:
            header = f.readline().rstrip('\n').split(',')
            gene_idx = {g: header.index(g) for g in TARGET if g in header}
            cluster_idx = header.index('cluster') if 'cluster' in header else 2
            for line in f:
                parts = line.rstrip('\n').split(',')
                if parts[cluster_idx] == 'beta':
                    ncell += 1
                    for g, gi in gene_idx.items():
                        gene_lists[g].append(float(parts[gi]))
        donor_beta_counts[tag] = gene_lists
        donor_beta_ncell[tag] = ncell
        print(f'{tag} ({DONOR_DISEASE[tag]}): beta cells = {ncell}')

    # donor-level aggregation (sum of counts)
    print('\n=== donor-level aggregation (sum of counts) ===')
    donor_sum = {g: {} for g in TARGET}
    for g in TARGET:
        row = []
        for tag in ['human1', 'human2', 'human3', 'human4']:
            s = sum(donor_beta_counts[tag][g])
            donor_sum[g][tag] = s
            row.append(f'{s:12.1f}')
        print(f'{g:10s} ' + ' '.join(row))

    out_csv = os.path.join(D, 'baron_donor_level_counts.csv')
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'human1_ND', 'human2_ND', 'human3_ND', 'human4_T2D'])
        for g in TARGET:
            w.writerow([g, donor_sum[g]['human1'], donor_sum[g]['human2'],
                        donor_sum[g]['human3'], donor_sum[g]['human4']])

    # directional log2FC (no statistical test)
    print('\n=== donor-level log2FC (3 ND donors vs 1 T2D donor, directional only) ===')
    out_fc = os.path.join(D, 'baron_donor_level_log2fc.csv')
    with open(out_fc, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'ND_mean_counts', 'T2D_counts', 'log2FC'])
        for g in TARGET:
            nd_vals = [donor_sum[g][t] for t in ['human1', 'human2', 'human3']]
            nd_mean = np.mean(nd_vals)
            t2d = donor_sum[g]['human4']
            lfc = np.log2((t2d + 0.5) / (nd_mean + 0.5))
            w.writerow([g, f'{nd_mean:.1f}', f'{t2d:.1f}', f'{lfc:.3f}'])
            print(f'{g:10s} {nd_mean:10.1f} {t2d:10.1f} {lfc:+8.3f}')

    print('T2D group has only 1 donor (human4); no statistical test is run, only a directional log2FC is reported.')


if __name__ == '__main__':
    main()
