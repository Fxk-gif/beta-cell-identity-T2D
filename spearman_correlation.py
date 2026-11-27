# -*- coding: utf-8 -*-
"""
Spearman rank-correlation robustness check for Reviewer 2 comment (3) -- real run (final)
================================================================
Data source: Supplementary Table S3 (table 3 of Supplementary_Materials_v2.docx)
        two columns of log2FC for the 15 common identity genes, kept to 3 decimals:
        GSE76895 (discovery cohort) / GSE76894 (validation cohort)
Standard: identical to the Pearson leave-one-gene-out analysis (same 15 genes, same two columns).

Key conclusions:
  Pearson  r   (all 15 genes) = 0.668    -> -0.19 after removing SLC2A2 (collapse)
  Spearman rho (all 15 genes) = 0.004    -> -0.23 after removing SLC2A2 (near zero/negative)
  => the rank correlation is essentially absent, further confirming that r=0.668 is driven entirely by the single high-leverage value SLC2A2,
     and does not support a "network-level coordinated downregulation".
"""
import math

# Table S3 exact values (3 decimals)
DATA = [
    ('SLC2A2', -1.418, -2.002), ('GLP1R', -0.409, -0.749), ('GCK', 0.062, -0.592),
    ('ABCC8', -0.061, -0.982), ('IAPP', -0.066, -0.788), ('NKX6-1', -0.452, -0.280),
    ('HNF1A', -0.320, -0.202), ('MAFB', -0.065, -0.262), ('HNF4A', -0.143, 0.136),
    ('NEUROD1', 0.058, -0.517), ('PAX6', 0.024, -0.485), ('SLC30A8', 0.005, -0.484),
    ('PCSK1', 0.114, -0.714), ('PCSK2', 0.190, -0.436), ('INS', -0.005, -0.254),
]
genes = [d[0] for d in DATA]
x = [d[1] for d in DATA]
y = [d[2] for d in DATA]
n = len(genes)


def pearson(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(len(a)))
    va = sum((a[i] - ma) ** 2 for i in range(len(a)))
    vb = sum((b[i] - mb) ** 2 for i in range(len(b)))
    return cov / math.sqrt(va * vb)


def rankdata(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[o[j + 1]] == v[o[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            r[o[k]] = avg
        i = j + 1
    return r


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def main():
    print('Data source: Supplementary Table S3 (15 genes, 3 decimals)')
    print('=' * 70)
    r_all = pearson(x, y)
    rho_all = spearman(x, y)
    print(f'all 15 genes   Pearson  r   = {r_all:.4f}')
    print(f'all 15 genes   Spearman rho = {rho_all:.4f}')
    print()
    idx = genes.index('SLC2A2')
    gx = [x[j] for j in range(n) if j != idx]
    gy = [y[j] for j in range(n) if j != idx]
    print(f'no SLC2A2     Pearson  r   = {pearson(gx, gy):.4f}')
    print(f'no SLC2A2     Spearman rho = {spearman(gx, gy):.4f}')
    print()
    others_r = [pearson([x[j] for j in range(n) if j != i], [y[j] for j in range(n) if j != i])
                for i in range(n) if genes[i] != 'SLC2A2']
    others_rho = [spearman([x[j] for j in range(n) if j != i], [y[j] for j in range(n) if j != i])
                  for i in range(n) if genes[i] != 'SLC2A2']
    print(f'leave-one-out Pearson  r   range: {min(others_r):.3f} ~ {max(others_r):.3f}')
    print(f'leave-one-out Spearman rho range: {min(others_rho):.3f} ~ {max(others_rho):.3f}')


if __name__ == '__main__':
    main()
