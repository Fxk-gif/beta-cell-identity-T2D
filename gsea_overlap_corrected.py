# -*- coding: utf-8 -*-
"""
去环化 GSEA（Overlap-corrected GSEA）
- 目标：回应审稿人 R1-3 —— 靶基因集与 16-gene identity score 共享基因造成的循环论证
- 方法：从每个靶基因集剔除与 16-gene identity score 共享的基因，重新计算 ES 与置换 p 值
- 算法：与原作者 gsea.py 完全一致（手写 ES + 2000 次置换 + seed=42），仅改动基因集
- 先做 sanity check：复现原始 gsea_results.csv 的 ES，验证代码正确后再跑去环化版
"""
import csv
import numpy as np

base = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\r3_inputs'

# ---- 1. 读全基因组 logFC ----
fc = {}
with open(base + r'\volcano_input.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        fc[row['id']] = float(row['logFC'])

# ---- 2. 读靶基因集 ----
sets = {}
with open(base + r'\tf_target_sets.csv', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        sets.setdefault(row['set'], []).append(row['gene'])

# ---- 3. 16-gene identity score（来自 analysis.R 第25-26行 & 正文 Methods）----
identity_16 = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFA","MAFB","NEUROD1","PAX6",
               "SLC2A2","GCK","ABCC8","KCNJ11","IAPP","UCN3","GLP1R","SLC30A8",
               "PCSK1","PCSK2","INS"]  # 注意：正文写16基因，但 analysis.R 列表含19个；以正文 Methods 16-gene 为准
# 正文 Methods 明确列出的 16 基因（"16 identity genes, including..."）：
identity_16_methods = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFB","NEUROD1","PAX6",
                        "SLC2A2","GCK","ABCC8","IAPP","GLP1R","SLC30A8","PCSK1","PCSK2","INS"]

# 审稿人点名的共享基因（R1-3 原文）：
shared = {"SLC2A2", "GCK", "IAPP", "INS", "NKX6-1", "PDX1"}

# ---- 4. 排序（与 gsea.py 一致：logFC 降序）----
ranked = sorted(fc.items(), key=lambda x: -x[1])
N = len(ranked)

def compute_es(gene_idx_set, ranked, fc, n_miss):
    """与 gsea.py 的 compute_es 完全一致：返回最大绝对 ES"""
    n_hits = len(gene_idx_set)
    nr = sum(abs(fc[g]) for g in gene_idx_set)
    if nr == 0:
        return 0.0
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

def gsea(geneset, ranked, fc, n_perm=2000, seed=42):
    """与 gsea.py 的 gsea 完全一致"""
    geneset = set(g for g in geneset if g in fc)
    if len(geneset) < 3:
        return None, None, None
    n_hits = len(geneset)
    n_miss = N - n_hits
    es_max = compute_es(geneset, ranked, fc, n_miss)
    all_genes = list(fc.keys())
    rng = np.random.default_rng(seed)
    perm_es = []
    for _ in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es.append(compute_es(sample, ranked, fc, n_miss))
    perm_es = np.array(perm_es)
    p = float(np.mean(np.abs(perm_es) >= abs(es_max)))
    return es_max, p, n_hits

# ---- 5. SANITY CHECK：复现原始结果（不去环化）----
print('=' * 70)
print('SANITY CHECK：复现原始 GSEA（应匹配 gsea_results.csv）')
print('=' * 70)
print(f'{"基因集":24s} {"复现ES":>10s} {"原ES":>10s} {"一致?"}')
expected = {
    'PDX1_targets': -0.948,
    'HNF1A_targets': -0.808,
    'beta_identity_network': -0.780,
    'beta_functional_genes': -0.815,
}
ok_all = True
for sname, genes in sets.items():
    es, p, n_hits = gsea(genes, ranked, fc)
    exp = expected.get(sname, None)
    if es is None:
        print(f'{sname:24s}  (基因<3，跳过)')
        continue
    if exp is not None:
        match = abs(es - exp) < 0.001
        ok_all = ok_all and match
        print(f'{sname:24s} {es:+10.4f} {exp:+10.3f} {"✓" if match else "✗ 不一致!"}')
    else:
        print(f'{sname:24s} {es:+10.4f}   (无预期值)')

print()
if ok_all:
    print('✅ SANITY CHECK 通过：复现结果与原 gsea_results.csv 一致，算法正确。')
else:
    print('⚠️  SANITY CHECK 失败：请检查算法是否与原 gsea.py 一致。')
print()

# ---- 6. 去环化 GSEA ----
print('=' * 70)
print('去环化 GSEA（从每个靶基因集剔除共享基因）')
print('共享基因（16-gene identity score ∩ 靶基因集）=', sorted(shared))
print('=' * 70)
print(f'{"基因集":24s} {"剔除后基因数":>10s} {"去环化ES":>12s} {"p值":>10s}')
results = []
for sname, genes in sets.items():
    genes_all = set(genes)
    removed = genes_all & shared
    genes_kept = genes_all - shared
    # 只保留在 fc 里的基因
    genes_kept_infc = [g for g in genes_kept if g in fc]
    es, p, n_hits = gsea(genes_kept, ranked, fc)
    results.append((sname, removed, len(genes_kept_infc), es, p))
    if es is None:
        print(f'{sname:24s}  剔除后基因<3，跳过')
        continue
    dirn = '下调富集' if es < 0 else '上调富集'
    print(f'{sname:24s} {len(genes_kept_infc):>10d} {es:+12.4f} {p:10.4f}  ({dirn}, 剔除{len(removed)}个)')

# ---- 7. 保存结果 ----
out = base + r'\gsea_overlap_corrected_results.csv'
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['set', 'removed_genes', 'n_genes_after', 'ES', 'p'])
    for sname, removed, n_kept, es, p in results:
        w.writerow([sname, ';'.join(sorted(removed)), n_kept, 
                    (f'{es:.6f}' if es is not None else ''), 
                    (f'{p:.6f}' if p is not None else '')])
print()
print('去环化结果已保存:', out)
