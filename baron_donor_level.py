# -*- coding: utf-8 -*-
"""
Baron (GSE84133) donor-level 重提取 —— 修复 cell-level 假重复问题

背景：审稿人1(R1-1)/审稿人2(R2-1)指出，Baron T2DM 组由单个供体组成，
cell-level Mann-Whitney 检验（FDR=2.4e-34 等）是假重复，样本量被严重高估。

本脚本：从 Baron 原始 counts 文件，保留 donor 归属（human1/2/3=ND, human4=T2D），
对每个 donor 的 β-cell 做 donor-level 聚合（每 donor 每基因的 counts 求和/均值），
输出 donor-level 的每基因 counts，供后续 donor-level 比较（T2D n=1，只报 log2FC）。

数据源：sc_data/GSM2230757_human1_umifm_counts.csv.gz ~ GSM2230760_human4_umifm_counts.csv.gz
（Baron et al. 2016, Cell Systems; human1-3 非糖尿病, human4 T2D）
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

# 身份基因（与正文一致，24 基因）
target = ['SLC2A2','IAPP','PDX1','NKX6-1','HNF1A','HNF4A','MAFA','MAFB','NEUROD1',
          'GCK','GLP1R','INS','ABCC8','KCNJ11','SLC30A8','PCSK2','CHL1',
          'GAD1','CALCA','TSHR','NHLH2','FOXE1','OPRD1','GRIN2D']

# 1) 提取每个 donor 的 β-cell 中每个目标基因的 counts 列表
donor_beta_counts = {}   # donor -> gene -> list[count]
donor_beta_ncell = {}    # donor -> β-cell 数
for tag, fn in files.items():
    path = os.path.join(D, fn)
    gene_lists = {g: [] for g in target}
    ncell = 0
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
        # 找 gene 列（第0列是 gene symbol? 需要确认列结构）
        # Baron 原始文件列结构: 通常 [cell, cluster, geneA, geneB, ...] 或 [gene, cell, ...]
        # 实测 baron_extract.py 里 parts[2]=='beta' 是 cluster，gene 用 header.index(g)
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
    print(f'{tag} ({donor_disease[tag]}): β-cell 数 = {ncell}')

# 2) donor-level 聚合：每 donor 每基因 counts 总和（pseudobulk sum）
print('\n=== donor-level 聚合 (counts 总和) ===')
print(f'{"gene":10s} ' + ' '.join(f'{d:>12s}' for d in ['human1(ND)','human2(ND)','human3(ND)','human4(T2D)']))
donor_sum = {g: {} for g in target}
for g in target:
    row = []
    for tag in ['human1','human2','human3','human4']:
        s = sum(donor_beta_counts[tag][g])
        donor_sum[g][tag] = s
        row.append(f'{s:12.1f}')
    print(f'{g:10s} ' + ' '.join(row))

# 3) 保存 donor-level 结果
out_csv = os.path.join(OUT, 'baron_donor_level_counts.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'human1_ND', 'human2_ND', 'human3_ND', 'human4_T2D'])
    for g in target:
        w.writerow([g, donor_sum[g]['human1'], donor_sum[g]['human2'],
                    donor_sum[g]['human3'], donor_sum[g]['human4']])

# 4) 计算 donor-level log2FC（ND 3 供体均值 vs T2D 1 供体），不报 p 值
print('\n=== donor-level log2FC (ND 3 供体 vs T2D 1 供体, 仅方向性) ===')
print(f'{"gene":10s} {"ND_mean":>10s} {"T2D":>10s} {"log2FC":>8s}')
out_fc = os.path.join(OUT, 'baron_donor_level_log2fc.csv')
with open(out_fc, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'ND_mean_counts', 'T2D_counts', 'log2FC'])
    for g in target:
        nd_vals = [donor_sum[g][t] for t in ['human1','human2','human3']]
        nd_mean = np.mean(nd_vals)
        t2d = donor_sum[g]['human4']
        # 用 +0.5 伪计数避免 log2(0)
        lfc = np.log2((t2d + 0.5) / (nd_mean + 0.5))
        w.writerow([g, f'{nd_mean:.1f}', f'{t2d:.1f}', f'{lfc:.3f}'])
        print(f'{g:10s} {nd_mean:10.1f} {t2d:10.1f} {lfc:+8.3f}')

print(f'\n结果已保存: {out_csv}, {out_fc}')
print('说明: T2D 组仅 1 供体(human4)，不做统计检验，只报方向性 log2FC（符合审稿人要求）')
