import gzip, csv
from collections import Counter

sc_dir = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\sc_data'

# 1. 读注释（GSM -> condition/subtype/donor），并按 GSM 顺序（即 Sample 顺序）
annot = {}
gsm_order = []
with open(sc_dir + r'\GSE81608_annot.csv', encoding='utf-8') as f:
    r = csv.DictReader(f)
    for row in r:
        gsm = row['GSM']
        gsm_order.append(gsm)
        annot[gsm] = row

print('注释样本数:', len(gsm_order))
print('condition 分布:', Counter(annot[g]['condition'] for g in gsm_order))
print('subtype 分布:', Counter(annot[g]['subtype'] for g in gsm_order))

# 2. 目标基因 Entrez ID
target = {
    'PDX1': 3651, 'NKX6-1': 4825, 'SLC2A2': 6514, 'GLP1R': 2740, 'HNF1A': 6927,
    'HNF4A': 3172, 'MAFA': 389692, 'MAFB': 9935, 'NEUROD1': 4760, 'PAX6': 5080,
    'GCK': 2645, 'ABCC8': 6833, 'KCNJ11': 3767, 'IAPP': 3375, 'SLC30A8': 169026,
    'INS': 3630, 'GAD1': 2571, 'GAD2': 2572, 'CALCA': 796, 'TSHR': 7253,
    'NHLH2': 4808, 'FOXE1': 2304, 'MYCN': 4613, 'CHL1': 10752, 'OPRD1': 4985,
    'GRIN2D': 2906, 'UCN3': 114131, 'PCSK1': 5122, 'PCSK2': 5126, 'NEUROG3': 50674,
}
id2sym = {str(v): k for k, v in target.items()}

# 3. 读取矩阵，提取目标基因表达（gene x sample）
# 矩阵：行=gene(Entrez), 列=Sample_1..1600
expr = {sym: None for sym in target}  # sym -> list of 1600 values
mat_file = sc_dir + r'\GSE81608_human_islets_rpkm.txt.gz'
with gzip.open(mat_file, 'rt') as f:
    header = f.readline().rstrip('\n').split('\t')
    sample_cols = header[1:]  # Sample_1 ...
    assert len(sample_cols) == 1600, f'样本列数 {len(sample_cols)}'
    for line in f:
        parts = line.rstrip('\n').split('\t')
        gid = parts[0]
        if gid in id2sym:
            sym = id2sym[gid]
            expr[sym] = [float(x) for x in parts[1:]]

missing = [s for s, v in expr.items() if v is None]
print('缺失基因:', missing)
print('成功提取基因数:', len([s for s, v in expr.items() if v is not None]))

# 4. 定位 beta 细胞样本索引
beta_idx = [i for i, g in enumerate(gsm_order) if annot[g]['subtype'] == 'beta']
print('beta 细胞数:', len(beta_idx))
# 检查 condition
beta_cond = [annot[gsm_order[i]]['condition'] for i in beta_idx]
print('beta condition:', Counter(beta_cond))

# 5. 提取 beta 细胞的身份基因表达矩阵
identity_genes = ['PDX1','NKX6-1','SLC2A2','GLP1R','HNF1A','HNF4A','MAFA','MAFB',
                  'NEUROD1','PAX6','GCK','ABCC8','KCNJ11','IAPP','SLC30A8','INS']
neuro_genes = ['GAD1','GAD2','CALCA','TSHR','NHLH2']

# 保存 beta 细胞表达子集
with open(sc_dir + r'\xin_beta_expr.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene'] + [f'S{i}' for i in beta_idx])
    for sym in list(target):
        if expr[sym] is not None:
            w.writerow([sym] + [expr[sym][i] for i in beta_idx])

# 保存 beta 注释
with open(sc_dir + r'\xin_beta_annot.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sample_idx', 'condition', 'donor'])
    for i in beta_idx:
        w.writerow([i, annot[gsm_order[i]]['condition'], annot[gsm_order[i]]['donor']])

print('\n已保存 xin_beta_expr.csv 和 xin_beta_annot.csv')
