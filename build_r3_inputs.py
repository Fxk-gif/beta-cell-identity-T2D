import csv, os

BASE = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\geo_data\processed'
R3IN = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44\r3_inputs'
os.makedirs(R3IN, exist_ok=True)

# ---- 1) 铁死亡基因文件 (KEGG hsa04216) ----
ferro = ['LPCAT3','ATG7','SAT2','CP','CYBB','ACSL1','ACSL3','ACSL4','ACSL6','SLC39A14',
         'SLC7A11','ALOX15','FTH1','FTL','GCLC','GCLM','GPX4','GSS','SLC40A1','HMOX1',
         'MAP1LC3BP1','MAP1LC3C','SLC11A2','PCBP1','PCBP2','ACSL5','STEAP3','PRNP','SAT1',
         'SLC39A8','MAP1LC3B2','SLC3A2','TF','TFRC','TP53','VDAC2','VDAC3','NCOA4',
         'MAP1LC3B','MAP1LC3A','FTMT','ATG5']
with open(os.path.join(R3IN, 'ferroptosis_genes.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gene'])
    for g in ferro:
        w.writerow([g])
print('ferroptosis_genes.csv:', len(ferro), 'genes')

# ---- 2) 训练集 GSE76895: 只保留 T2D/ND ----
grp = {}
with open(os.path.join(BASE, 'GSE76895_groups_full.csv'), encoding='utf-8') as f:
    for row in csv.DictReader(f):
        grp[row['sample']] = row['group']

keep = [s for s, g in grp.items() if g in ('T2D', 'ND')]
src = os.path.join(BASE, 'GSE76895_expr_gene.csv')
dst = os.path.join(R3IN, 'GSE76895_expr.csv')
with open(src, encoding='utf-8') as f, open(dst, 'w', newline='', encoding='utf-8') as out:
    rd = csv.reader(f)
    header = next(rd)
    idx = {s: i for i, s in enumerate(header)}
    cols = [0] + [idx[s] for s in keep]  # 基因列 + 保留样本
    w = csv.writer(out)
    w.writerow([header[0]] + keep)
    n = 0
    for row in rd:
        w.writerow([row[i] for i in cols])
        n += 1
print('GSE76895_expr.csv:', n, 'genes x', len(keep), 'samples (T2D+ND)')

with open(os.path.join(R3IN, 'GSE76895_groups.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['sample', 'group'])
    for s in keep:
        w.writerow([s, grp[s]])
print('GSE76895_groups.csv: written')

# ---- 3) 验证集 GSE76894: 全部103例 ----
grp94 = {}
with open(os.path.join(BASE, 'GSE76894_groups_full.csv'), encoding='utf-8') as f:
    for row in csv.DictReader(f):
        grp94[row['sample']] = row['group']
src = os.path.join(BASE, 'GSE76894_expr_gene.csv')
dst = os.path.join(R3IN, 'GSE76894_expr.csv')
with open(src, encoding='utf-8') as f, open(dst, 'w', newline='', encoding='utf-8') as out:
    rd = csv.reader(f)
    header = next(rd)
    samples = header[1:]
    w = csv.writer(out)
    w.writerow(header)
    n = 0
    for row in rd:
        w.writerow(row)
        n += 1
print('GSE76894_expr.csv:', n, 'genes x', len(samples), 'samples')

with open(os.path.join(R3IN, 'GSE76894_groups.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['sample', 'group'])
    for s in samples:
        w.writerow([s, grp94[s]])
print('GSE76894_groups.csv: written')

# ---- 校验 ----
t2d = sum(1 for s in keep if grp[s]=='T2D'); nd = sum(1 for s in keep if grp[s]=='ND')
t2d94 = sum(1 for g in grp94.values() if g=='T2D'); nd94 = sum(1 for g in grp94.values() if g=='ND')
print(f'TRAIN GSE76895: T2D={t2d} ND={nd}')
print(f'VALID GSE76894: T2D={t2d94} ND={nd94}')
print('ALL DONE')
