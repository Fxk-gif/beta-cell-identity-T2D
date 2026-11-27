import gzip, os, csv
import paths
import numpy as np
from scipy import stats
from collections import Counter

d = paths.sc_data()

# donor -> disease status
donor_disease = {
    'human1': 'ND', 'human2': 'ND', 'human3': 'ND', 'human4': 'T2D'
}
files = {
    'human1': 'GSM2230757_human1_umifm_counts.csv.gz',
    'human2': 'GSM2230758_human2_umifm_counts.csv.gz',
    'human3': 'GSM2230759_human3_umifm_counts.csv.gz',
    'human4': 'GSM2230760_human4_umifm_counts.csv.gz',
}

# target genes
target = ['SLC2A2','IAPP','PDX1','NKX6-1','HNF1A','HNF4A','MAFA','MAFB','NEUROD1',
          'GCK','GLP1R','INS','ABCC8','KCNJ11','SLC30A8','PCSK2','CHL1',
          'GAD1','CALCA','TSHR','NHLH2','FOXE1','OPRD1','GRIN2D']

# extract target-gene counts from each donor's beta cells
all_beta = {}  # gene -> list of (donor, count)
for tag, fn in files.items():
    path = os.path.join(d, fn)
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
        gene_idx = {g: header.index(g) for g in target if g in header}
        # cluster column index = 2
        for line in f:
            parts = line.rstrip('\n').split(',')
            if parts[2] == 'beta':
                for g, gi in gene_idx.items():
                    all_beta.setdefault(g, []).append((tag, int(parts[gi])))

# convert to numpy, grouped by disease
gene_expr = {}  # gene -> {ND: array, T2D: array}
for g in target:
    rows = all_beta.get(g, [])
    nd = np.array([v for (d_, v) in rows if donor_disease[d_] == 'ND'], dtype=float)
    t2d = np.array([v for (d_, v) in rows if donor_disease[d_] == 'T2D'], dtype=float)
    gene_expr[g] = {'ND': nd, 'T2D': t2d}

# save the beta expression as CSV (for R plotting)
out = os.path.join(d, 'baron_beta_expr.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene'] + [f'ND_{i}' for i in range(2030)] + [f'T2D_{i}' for i in range(495)])
    for g in target:
        nd = gene_expr[g]['ND']
        t2d = gene_expr[g]['T2D']
        # pad lengths (the beta count per gene should match, all from the same beta cells)
        w.writerow([g] + list(nd) + list(t2d))

# statistical test
print(f'{"gene":10s} {"ND_mean":>9s} {"T2D_mean":>9s} {"log2FC":>8s} {"p":>10s}')
results = []
for g in target:
    nd = gene_expr[g]['ND']
    t2d = gene_expr[g]['T2D']
    if len(nd) == 0 or len(t2d) == 0:
        continue
    ndm, t2dm = nd.mean(), t2d.mean()
    # Mann-Whitney U test
    p = stats.mannwhitneyu(t2d, nd, alternative='two-sided').pvalue
    lfc = np.log2((t2dm + 0.1) / (ndm + 0.1))
    results.append((g, ndm, t2dm, lfc, p, len(nd), len(t2d)))
    print(f'{g:10s} {ndm:9.2f} {t2dm:9.2f} {lfc:+8.3f} {p:10.3e}')

# FDR correction
from statsmodels.stats.multitest import multipletests
ps = [r[4] for r in results]
rej, fdr, _, _ = multipletests(ps, method='fdr_bh')
print('\n=== after FDR correction ===')
for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
    flag = '**significant**' if f < 0.05 else ''
    print(f'{g:10s} log2FC={lfc:+7.3f} p={p:.3e} FDR={f:.3e} {flag}')

# save the statistics
with open(os.path.join(d, 'baron_beta_stats.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'ND_mean', 'T2D_mean', 'log2FC', 'p', 'FDR', 'ND_n', 'T2D_n'])
    for (g, ndm, t2dm, lfc, p, n1, n2), fdr_v in zip(results, fdr):
        w.writerow([g, f'{ndm:.3f}', f'{t2dm:.3f}', f'{lfc:.3f}', f'{p:.3e}', f'{fdr_v:.3e}', n1, n2])

print('\nResults saved: baron_beta_expr.csv, baron_beta_stats.csv')
n_nd = sum(len(gene_expr[g]["ND"]) for g in target) // len(target)
n_t2d = sum(len(gene_expr[g]["T2D"]) for g in target) // len(target)
print('beta cells: ND=%d T2D=%d' % (n_nd, n_t2d))
