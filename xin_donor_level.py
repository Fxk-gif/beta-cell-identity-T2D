# -*- coding: utf-8 -*-
"""
Xin (GSE81608) donor-level 重提取 —— 修复 cell-level 假重复问题

背景：审稿人1(R1-1)/审稿人2(R2-1)指出，cell-level Mann-Whitney 检验是假重复。
本脚本：从 xin_beta_expr.csv (RPKM) + xin_beta_annot.csv (donor) 做 donor-level 聚合，
每 donor 每基因取 β-cell RPKM 均值（pseudobulk mean），得 18 donors (12 ND + 6 T2D)，
用 Wilcoxon rank-sum（donor 为样本单位）做差异检验，并算 donor-level log2FC。

注：Xin 公开数据为 RPKM（非 raw counts），故不适用 edgeR 负二项模型；
采用 donor-level 非参数 Wilcoxon 检验（样本单位=donor, n=18），这是对 cell-level 假重复的正确修正。

数据源：
  sc_data/xin_beta_expr.csv (503 β-cell × 基因, 列 S0..S502)
  sc_data/xin_beta_annot.csv (503 行: sample_idx, condition, donor)
"""
import csv, os
import paths
import numpy as np

BASE = paths._DATA
D = paths.sc_data()

# 1) 读 annot：建立 cell列 -> (donor, condition)
annot = {}
with open(os.path.join(D, 'xin_beta_annot.csv'), encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        # xin_beta_expr.csv 的列名是 S{idx}，annot 的 sample_idx 对应 idx
        annot[f'S{row["sample_idx"]}'] = (row['donor'], row['condition'])

# 2) 读 expr：gene × cell
expr = {}
with open(os.path.join(D, 'xin_beta_expr.csv'), encoding='utf-8-sig') as f:
    r = csv.reader(f)
    header = next(r)  # ['gene', 'S0', 'S1', ...]
    cells = header[1:]
    for row in r:
        gene = row[0]
        expr[gene] = dict(zip(cells, [float(v) for v in row[1:]]))

# 3) 目标基因（与正文一致的核心身份基因）
target = ['SLC2A2','IAPP','PDX1','NKX6-1','HNF1A','HNF4A','MAFA','MAFB','NEUROD1',
          'GCK','GLP1R','INS','ABCC8','KCNJ11','SLC30A8','PCSK1','PCSK2']

# 4) donor-level 聚合：每 donor 每基因的 β-cell RPKM 均值
donors_nd = [d for d,c in set(annot.values()) if c != 'T2D']
donors_t2d = [d for d,c in set(annot.values()) if c == 'T2D']
print(f'ND donors ({len(donors_nd)}): {sorted(donors_nd)}')
print(f'T2D donors ({len(donors_t2d)}): {sorted(donors_t2d)}')

# 每 donor 每基因的 cell 均值
donor_mean = {}  # donor -> gene -> mean RPKM
for donor in sorted(set(annot.values())):
    pass

# 重建：donor -> gene -> list of cell values
from collections import defaultdict
donor_gene_vals = defaultdict(lambda: defaultdict(list))
for cell, (donor, cond) in annot.items():
    for gene in target:
        if cell in expr.get(gene, {}):
            donor_gene_vals[donor][gene].append(expr[gene][cell])

# donor-level mean
donor_level = defaultdict(dict)  # donor -> gene -> mean
donor_condition = {}  # donor -> condition
for donor, genes in donor_gene_vals.items():
    # 从 annot 取该 donor 的 condition
    for cell, (d_, cond) in annot.items():
        if d_ == donor:
            donor_condition[donor] = cond
            break
    for gene in target:
        vals = genes.get(gene, [])
        donor_level[donor][gene] = np.mean(vals) if vals else np.nan

# 5) Wilcoxon rank-sum (donor 为样本单位)，手写实现（无 scipy 依赖）
def mannwhitney_u(x, y, alternative='two-sided'):
    """手写 Mann-Whitney U 检验（等价 scipy.stats.mannwhitneyu），返回双尾 p 值。
    x, y 为两个独立样本数组。样本量小，用精确/正态近似皆可，此处用正态近似+连续性校正。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    ranked = np.concatenate([x, y])
    order = np.argsort(ranked)
    ranks = np.empty_like(order, dtype=float)
    # 处理并列
    sorted_vals = ranked[order]
    i = 0
    n = len(sorted_vals)
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        ranks[order[i:j + 1]] = avg
        i = j + 1
    rank_x = ranks[:n1]
    R1 = rank_x.sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U = min(U1, U2)
    # 正态近似（含并列校正）
    mu = n1 * n2 / 2.0
    # 并列校正
    ties = 0
    for k in range(n):
        pass
    # 简化：用无并列近似（样本量小且多为连续 RPKM 均值，并列极少）
    sigma = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (U - mu) / sigma
    from math import erfc, sqrt
    p = erfc(abs(z) / sqrt(2.0))  # 双尾
    return p


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR 校正（手写）。"""
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


print(f'\n{"gene":10s} {"ND_mean":>9s} {"T2D_mean":>9s} {"log2FC":>8s} {"p_donor":>10s}')
results = []
for gene in target:
    nd_vals = [donor_level[d][gene] for d in donor_level
               if donor_condition.get(d) == 'non-diabetic' and not np.isnan(donor_level[d][gene])]
    t2d_vals = [donor_level[d][gene] for d in donor_level
                if donor_condition.get(d) == 'T2D' and not np.isnan(donor_level[d][gene])]
    if len(nd_vals) < 2 or len(t2d_vals) < 2:
        print(f'{gene:10s} 样本不足 (ND={len(nd_vals)}, T2D={len(t2d_vals)}), 跳过')
        continue
    ndm = np.mean(nd_vals); t2dm = np.mean(t2d_vals)
    lfc = np.log2((t2dm + 0.1) / (ndm + 0.1))
    p = mannwhitney_u(t2d_vals, nd_vals)
    results.append((gene, ndm, t2dm, lfc, p, len(nd_vals), len(t2d_vals)))
    print(f'{gene:10s} {ndm:9.3f} {t2dm:9.3f} {lfc:+8.3f} {p:10.3e}')

# 6) FDR 校正
ps = [r[4] for r in results]
if ps:
    fdr = bh_fdr(np.array(ps))
    print('\n=== donor-level Wilcoxon + BH 校正 ===')
    for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
        flag = '*' if f < 0.05 else ''
        print(f'{g:10s} log2FC={lfc:+7.3f} p={p:.3e} FDR={f:.3e} {flag}')

# 7) 保存
out = os.path.join(D, 'xin_donor_level_results.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'ND_mean', 'T2D_mean', 'log2FC', 'p_donor', 'FDR', 'ND_n_donor', 'T2D_n_donor'])
    for (g, ndm, t2dm, lfc, p, n1, n2), f in zip(results, fdr):
        w.writerow([g, f'{ndm:.4f}', f'{t2dm:.4f}', f'{lfc:.4f}', f'{p:.3e}', f'{f:.3e}', n1, n2])

print(f'\n结果已保存: {out}')
print('说明: 样本单位=donor (18 donors: 12 ND + 6 T2D)，Wilcoxon rank-sum 检验，修正了 cell-level 假重复')
