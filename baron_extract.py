import gzip, os, csv
import numpy as np
from scipy import stats
from collections import Counter

d = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\sc_data'

# 供体 -> 疾病状态
donor_disease = {
    'human1': 'ND', 'human2': 'ND', 'human3': 'ND', 'human4': 'T2D'
}
files = {
    'human1': 'GSM2230757_human1_umifm_counts.csv.gz',
    'human2': 'GSM2230758_human2_umifm_counts.csv.gz',
    'human3': 'GSM2230759_human3_umifm_counts.csv.gz',
    'human4': 'GSM2230760_human4_umifm_counts.csv.gz',
}

# 目标基因
target = ['SLC2A2','IAPP','PDX1','NKX6-1','HNF1A','HNF4A','MAFA','MAFB','NEUROD1',
          'GCK','GLP1R','INS','ABCC8','KCNJ11','SLC30A8','PCSK2','CHL1',
          'GAD1','CALCA','TSHR','NHLH2','FOXE1','OPRD1','GRIN2D']

# 提取每个供体 beta 细胞的目标基因 counts
all_beta = {}  # gene -> list of (donor, count)
for tag, fn in files.items():
    path = os.path.join(d, fn)
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
        gene_idx = {g: header.index(g) for g in target if g in header}
        # cluster 列索引 = 2
        for line in f:
            parts = line.rstrip('\n').split(',')
            if parts[2] == 'beta':
                for g, gi in gene_idx.items():
                    all_beta.setdefault(g, []).append((tag, int(parts[gi])))

# 转成 numpy，按疾病分组
gene_expr = {}  # gene -> {ND: array, T2D: array}
for g in target:
    rows = all_beta.get(g, [])
    nd = np.array([v for (d_, v) in rows if donor_disease[d_] == 'ND'], dtype=float)
    t2d = np.array([v for (d_, v) in rows if donor_disease[d_] == 'T2D'], dtype=float)
    gene_expr[g] = {'ND': nd, 'T2D': t2d}

# 保存 beta 表达为 CSV（供 R 绘图）
out = os.path.join(d, 'baron_beta_expr.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene'] + [f'ND_{i}' for i in range(2030)] + [f'T2D_{i}' for i in range(495)])
    for g in target:
        nd = gene_expr[g]['ND']
        t2d = gene_expr[g]['T2D']
        # 补齐长度（不同基因的beta数应该一样，因为都来自同一批beta细胞）
        w.writerow([g] + list(nd) + list(t2d))

# 统计检验
print(f'{"gene":10s} {"ND_mean":>9s} {"T2D_mean":>9s} {"log2FC":>8s} {"p":>10s}')
results = []
for g in target:
    nd = gene_expr[g]['ND']
    t2d = gene_expr[g]['T2D']
    if len(nd) == 0 or len(t2d) == 0:
        continue
    ndm, t2dm = nd.mean(), t2d.mean()
    # Mann-Whitney U 检验
    p = stats.mannwhitneyu(t2d, nd, alternative='two-sided').pvalue
    lfc = np.log2((t2dm + 0.1) / (ndm + 0.1))
    results.append((g, ndm, t2dm, lfc, p, len(nd), len(t2d)))
    print(f'{g:10s} {ndm:9.2f} {t2dm:9.2f} {lfc:+8.3f} {p:10.3e}')

# FDR 校正
from statsmodels.stats.multitest import multipletests
ps = [r[4] for r in results]
rej, fdr, _, _ = multipletests(ps, method='fdr_bh')
print('\n=== FDR 校正后 ===')
for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
    flag = '**显著**' if f < 0.05 else ''
    print(f'{g:10s} log2FC={lfc:+7.3f} p={p:.3e} FDR={f:.3e} {flag}')

# 保存统计结果
with open(os.path.join(d, 'baron_beta_stats.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'ND_mean', 'T2D_mean', 'log2FC', 'p', 'FDR', 'ND_n', 'T2D_n'])
    for (g, ndm, t2dm, lfc, p, n1, n2), fdr_v in zip(results, fdr):
        w.writerow([g, f'{ndm:.3f}', f'{t2dm:.3f}', f'{lfc:.3f}', f'{p:.3e}', f'{fdr_v:.3e}', n1, n2])

print('\n结果已保存: baron_beta_expr.csv, baron_beta_stats.csv')
n_nd = sum(len(gene_expr[g]["ND"]) for g in target) // len(target)
n_t2d = sum(len(gene_expr[g]["T2D"]) for g in target) // len(target)
print('beta细胞数: ND=%d T2D=%d' % (n_nd, n_t2d))
