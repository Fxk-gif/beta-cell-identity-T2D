# -*- coding: utf-8 -*-
"""
审稿人2 意见(3) 的 Spearman 秩相关稳健性检验 —— 真实运行（最终版）
================================================================
数据源：Supplementary Table S3（Supplementary_Materials_v2.docx 表格3）
        15 个共同身份基因的两列 log₂FC，保留 3 位小数：
        GSE76895 (发现队列) / GSE76894 (验证队列)
口径：与 Pearson leave-one-gene-out 分析完全一致（同一 15 基因、同一两列）。

关键结论：
  Pearson  r   (全15基因) = 0.8283   -> 去掉 SLC2A2 后 0.0761（崩塌）
  Spearman rho (全15基因) = 0.1055   -> 去掉 SLC2A2 后 -0.1013（接近零/负）
  => 秩相关完全不相关，进一步证实 r=0.828 完全由 SLC2A2 一个高杠杆值驱动，
     不支持"网络级协调下调"。
"""
import math

# Table S3 精确值（3 位小数）
DATA = [
    ('SLC2A2', -0.252, -0.474), ('GLP1R', -0.072, -0.143), ('GCK', 0.011, -0.120),
    ('ABCC8', -0.006, -0.120), ('IAPP', -0.007, -0.090), ('NKX6-1', -0.096, -0.065),
    ('HNF1A', -0.059, -0.038), ('MAFB', -0.007, -0.029), ('HNF4A', -0.019, 0.017),
    ('NEUROD1', 0.007, -0.070), ('PAX6', 0.003, -0.060), ('SLC30A8', 0.001, -0.055),
    ('PCSK1', 0.013, -0.087), ('PCSK2', 0.022, -0.056), ('INS', -0.001, -0.027),
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
    print('数据源: Supplementary Table S3 (15 基因, 3 位小数)')
    print('=' * 70)
    r_all = pearson(x, y)
    rho_all = spearman(x, y)
    print(f'全 15 基因   Pearson  r   = {r_all:.4f}')
    print(f'全 15 基因   Spearman rho = {rho_all:.4f}')
    print()
    idx = genes.index('SLC2A2')
    gx = [x[j] for j in range(n) if j != idx]
    gy = [y[j] for j in range(n) if j != idx]
    print(f'去 SLC2A2    Pearson  r   = {pearson(gx, gy):.4f}')
    print(f'去 SLC2A2    Spearman rho = {spearman(gx, gy):.4f}')
    print()
    others_r = [pearson([x[j] for j in range(n) if j != i], [y[j] for j in range(n) if j != i])
                for i in range(n) if genes[i] != 'SLC2A2']
    others_rho = [spearman([x[j] for j in range(n) if j != i], [y[j] for j in range(n) if j != i])
                  for i in range(n) if genes[i] != 'SLC2A2']
    print(f'去其他基因   Pearson  r   范围: {min(others_r):.3f} ~ {max(others_r):.3f}')
    print(f'去其他基因   Spearman rho 范围: {min(others_rho):.3f} ~ {max(others_rho):.3f}')


if __name__ == '__main__':
    main()
