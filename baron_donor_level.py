# -*- coding: utf-8 -*-
"""
Baron (GSE84133) donor-level re-extraction -- fixes the cell-level pseudoreplication issue

Background: Reviewers R1-1/R2-1 noted that the Baron T2DM group consists of a single donor, and the
cell-level Mann-Whitney test (FDR=2.4e-34 etc.) is pseudoreplicated with a greatly inflated sample size.

This script: reads the Baron raw counts, keeps the donor assignment (human1/2/3=ND, human4=T2D),
aggregates beta-cells per donor (sum/mean of counts per donor per gene),
and outputs donor-level per-gene counts for downstream donor-level comparison (T2D n=1, so only log2FC is reported).

Data sources: sc_data/GSM2230757_human1_umifm_counts.csv.gz ~ GSM2230760_human4_umifm_counts.csv.gz
(Baron et al. 2016, Cell Systems; human1-3 non-diabetic, human4 T2D)
"""
import gzip, os, csv
import paths
import numpy as np

BASE = paths._DATA
D = paths.sc_data()
OUT = paths.sc_data()

donor_disease = {'human1': 'ND', 'human2': 'ND', 'human3': 'ND', 'human4': 'T2D'}
files = {
    'human1': 'GSM2230757_human1_umifm_counts.csv.gz',
    'human2': 'GSM2230758_human2_umifm_counts.csv.gz',
    'human3': 'GSM2230759_human3_umifm_counts.csv.gz',
    'human4': 'GSM2230760_human4_umifm_counts.csv.gz',
}

# identity genes (24 genes, matching the manuscript)
target = ['SLC2A2','IAPP','PDX1','NKX6-1','HNF1A','HNF4A','MAFA','MAFB','NEUROD1',
          'GCK','GLP1R','INS','ABCC8','KCNJ11','SLC30A8','PCSK2','CHL1',
          'GAD1','CALCA','TSHR','NHLH2','FOXE1','OPRD1','GRIN2D']

# 1) extract, per donor, the count list of each target gene across beta-cells
donor_beta_counts = {}   # donor -> gene -> list[count]
donor_beta_ncell = {}    # donor -> number of beta-cells
for tag, fn in files.items():
    path = os.path.join(D, fn)
    gene_lists = {g: [] for g in target}
    ncell = 0
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
        # locate the gene column (column 0 is the gene symbol? verify the column layout)
        # Baron raw file layout: usually [cell, cluster, geneA, geneB, ...] or [gene, cell, ...]
        # in baron_extract.py, parts[2]=='beta' is the cluster; genes use header.index(g)
        gene_idx = {g: header.index(g) for g in target if g in header}
        cluster_idx = header.index('cluster') if 'cluster' in header else 2
        for line in f:
            parts = line.rstrip('\n').split(',')
            if parts[cluster_idx] == 'beta':
                ncell += 1
                for g, gi in gene_idx.items():
                    gene_lists[g].append(float(parts[gi]))
    donor_beta_counts[tag] = gene_lists
    donor_beta_ncell[tag] = ncell
    print(f'{tag} ({donor_disease[tag]}): beta-cells = {ncell}')

# 2) donor-level aggregation: sum of counts per donor per gene (pseudobulk sum)
print('\n=== donor-level aggregation (sum of counts) ===')
print(f'{"gene":10s} ' + ' '.join(f'{d:>12s}' for d in ['human1(ND)','human2(ND)','human3(ND)','human4(T2D)']))
donor_sum = {g: {} for g in target}
for g in target:
    row = []
    for tag in ['human1','human2','human3','human4']:
        s = sum(donor_beta_counts[tag][g])
        donor_sum[g][tag] = s
        row.append(f'{s:12.1f}')
    print(f'{g:10s} ' + ' '.join(row))

# 3) save donor-level results
out_csv = os.path.join(OUT, 'baron_donor_level_counts.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'human1_ND', 'human2_ND', 'human3_ND', 'human4_T2D'])
    for g in target:
        w.writerow([g, donor_sum[g]['human1'], donor_sum[g]['human2'],
                    donor_sum[g]['human3'], donor_sum[g]['human4']])

# 4) donor-level log2FC (3 ND donors' mean vs 1 T2D donor), no p-value reported
print('\n=== donor-level log2FC (3 ND donors vs 1 T2D donor, directional only) ===')
print(f'{"gene":10s} {"ND_mean":>10s} {"T2D":>10s} {"log2FC":>8s}')
out_fc = os.path.join(OUT, 'baron_donor_level_log2fc.csv')
with open(out_fc, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'ND_mean_counts', 'T2D_counts', 'log2FC'])
    for g in target:
        nd_vals = [donor_sum[g][t] for t in ['human1','human2','human3']]
        nd_mean = np.mean(nd_vals)
        t2d = donor_sum[g]['human4']
        # use a +0.5 pseudocount to avoid log2(0)
        lfc = np.log2((t2d + 0.5) / (nd_mean + 0.5))
        w.writerow([g, f'{nd_mean:.1f}', f'{t2d:.1f}', f'{lfc:.3f}'])
        print(f'{g:10s} {nd_mean:10.1f} {t2d:10.1f} {lfc:+8.3f}')

print(f'\nResults saved: {out_csv}, {out_fc}')
print('Note: the T2D group has a single donor (human4); no statistical test, only a directional log2FC (as the reviewer requested)')
