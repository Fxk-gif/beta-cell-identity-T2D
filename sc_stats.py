import csv
import numpy as np
from scipy import stats

sc_dir = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\sc_data'

# 读 beta 表达
expr = {}
gene_order = []
with open(sc_dir + r'\xin_beta_expr.csv', encoding='utf-8') as f:
    r = csv.reader(f)
    header = next(r)
    n_cells = len(header) - 1
    for row in r:
        sym = row[0]
        gene_order.append(sym)
        expr[sym] = np.array([float(x) for x in row[1:]])

# 读注释
cond = []
with open(sc_dir + r'\xin_beta_annot.csv', encoding='utf-8') as f:
    r = csv.DictReader(f)
    for row in r:
        cond.append(row['condition'])
cond = np.array(cond)  # 'T2D' or 'non-diabetic'

print('beta 细胞数:', len(cond), 'T2D:', sum(cond=='T2D'), 'ND:', sum(cond=='non-diabetic'))

# 分组索引
t2d = cond == 'T2D'
nd = cond == 'non-diabetic'

# 身份基因与神经内分泌基因
identity_genes = ['PDX1','NKX6-1','SLC2A2','GLP1R','HNF1A','HNF4A','MAFA','MAFB',
                  'NEUROD1','PAX6','GCK','ABCC8','KCNJ11','IAPP','SLC30A8','INS']
neuro_genes = ['GAD1','GAD2','CALCA','TSHR','NHLH2']

# 1. 身份评分（z-score 平均，跨基因）
def zscore(x):
    return (x - np.mean(x)) / np.std(x)

id_mat = np.array([zscore(expr[g]) for g in identity_genes])
identity_score = id_mat.mean(axis=0)

# 2. 统计检验
wi = stats.mannwhitneyu(identity_score[t2d], identity_score[nd], alternative='two-sided')
print('\n=== β身份评分（单细胞 Xin）===')
print(f'Mann-Whitney U p = {wi.pvalue:.3e}')
print(f'ND 均值: {identity_score[nd].mean():.3f}, T2D 均值: {identity_score[t2d].mean():.3f}')

# 3. 单基因检验
print('\n=== 单基因 T2D vs ND（β细胞）===')
print(f'{"gene":10s} {"ND_mean":>9s} {"T2D_mean":>9s} {"log2FC":>8s} {"p":>12s} {"padj":>12s}')
results = []
for g in gene_order:
    nd_vals = expr[g][nd]
    t2d_vals = expr[g][t2d]
    p = stats.mannwhitneyu(t2d_vals, nd_vals, alternative='two-sided').pvalue
    nd_m = nd_vals.mean()
    t2d_m = t2d_vals.mean()
    lfc = np.log2((t2d_m + 0.1) / (nd_m + 0.1))
    results.append((g, nd_m, t2d_m, lfc, p))

from scipy.stats import false_discovery_control
ps = np.array([r[4] for r in results])
padjs = false_discovery_control(ps, method='bh')
for (g, nd_m, t2d_m, lfc, p), padj in zip(results, padjs):
    print(f'{g:10s} {nd_m:9.2f} {t2d_m:9.2f} {lfc:+8.3f} {p:12.2e} {padj:12.2e}')

# 4. 保存完整结果
with open(sc_dir + r'\xin_beta_results.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene','ND_mean','T2D_mean','log2FC','p','padj'])
    for (g, nd_m, t2d_m, lfc, p), padj in zip(results, padjs):
        w.writerow([g, nd_m, t2d_m, lfc, p, padj])

# 保存身份评分
with open(sc_dir + r'\xin_identity_score.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['condition','identity_score'])
    for i in range(len(cond)):
        w.writerow([cond[i], identity_score[i]])

print('\n结果已保存')
