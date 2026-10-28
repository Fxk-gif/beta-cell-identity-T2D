import csv, os
import numpy as np
from scipy import stats

base = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\r3_inputs'

# 读全基因 logFC
fc = {}
with open(base + r'\volcano_input.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fc[row['id']] = float(row['logFC'])

# 读靶基因集
sets = {}
with open(base + r'\tf_target_sets.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sets.setdefault(row['set'], []).append(row['gene'])

# 按 logFC 降序排序（GSEA 标准：正 logFC 在前）
ranked = sorted(fc.items(), key=lambda x: -x[1])  # (gene, logFC) 降序
N = len(ranked)

def compute_es(gene_idx_set, ranked, fc, n_miss):
    """给定命中基因的index集合，计算最大ES（向量化）"""
    n_hits = len(gene_idx_set)
    nr = sum(abs(fc[g]) for g in gene_idx_set)
    es = 0.0
    es_max = 0.0
    for g, lfc in ranked:
        if g in gene_idx_set:
            es += abs(lfc) / nr
        else:
            es -= 1.0 / n_miss
        if abs(es) > abs(es_max):
            es_max = es
    return es_max

def gsea(geneset, ranked, n_perm=2000):
    """手写 GSEA，返回 ES 和 p 值"""
    geneset = set(g for g in geneset if g in fc)
    if len(geneset) < 3:
        return None, None, None
    n_hits = len(geneset)
    n_miss = N - n_hits
    # 真实 ES
    es_max = compute_es(geneset, ranked, fc, n_miss)
    # 置换检验
    all_genes = list(fc.keys())
    perm_es = []
    rng = np.random.default_rng(42)
    for _ in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es.append(compute_es(sample, ranked, fc, n_miss))
    perm_es = np.array(perm_es)
    p = float(np.mean(np.abs(perm_es) >= abs(es_max)))
    if es_max >= 0:
        p_up = float(np.mean(perm_es >= es_max))
    else:
        p_up = float(np.mean(perm_es <= es_max))
    return es_max, p, n_hits

print('=== GSEA 富集检验（全基因组 logFC 排序）===')
print(f'{"基因集":22s} {"ES":>8s} {"p(双侧)":>10s} {"n_hits":>7s}  方向')
results = []
for sname, genes in sets.items():
    es, p, n_hits = gsea(genes, ranked)
    if es is None:
        continue
    direction = '下调富集' if es < 0 else '上调富集'
    results.append((sname, es, p, n_hits))
    print(f'{sname:22s} {es:+8.3f} {p:10.4f} {n_hits:7d}  {direction}')

# 保存结果
with open(base + r'\gsea_results.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['set', 'ES', 'p', 'n_hits'])
    for r in results:
        w.writerow(list(r))
print('\nGSEA 结果已保存: gsea_results.csv')
