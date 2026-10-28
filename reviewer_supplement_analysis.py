# -*- coding: utf-8 -*-
"""
审稿意见补充分析：
#1 单细胞"精简协调评分"——用真正胞内下调的基因（SLC2A2 + IAPP + INS）构建精简身份评分，
   检验其在单细胞层面是否显著（回应"16基因评分 p=0.79 不显著 vs 协调侵蚀"的张力）。
#2 progenitor/disallowed 基因评分（bulk）——验证 T2DM 中反向更低（支撑"成熟耗竭而非去分化"）。
所有统计口径与正文一致：bulk 用 Wilcoxon（scipy mannwhitneyu，与 R wilcox.test 渐近一致），
单细胞用 Mann-Whitney U。z-score 平均法跨基因。
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44'
SC = BASE + r'\sc_data'
R3 = BASE + r'\r3_inputs'
OUT = BASE + r'\figures'
import os
os.makedirs(OUT, exist_ok=True)

def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)  # 与 R sd() 一致（n-1）
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)

def mw_p(a, b):
    return stats.mannwhitneyu(a, b, alternative='two-sided').pvalue

results = {}

# ============================================================
# 分析 #1：单细胞精简协调评分（SLC2A2 + IAPP + INS）
# ============================================================
core_cell_genes = ['SLC2A2', 'IAPP', 'INS']  # 真正胞内下调的功能基因

for name, expr_f, annot_f, cond_col, t2d_vals, nd_vals in [
    ('Xin', 'xin_beta_expr.csv', 'xin_beta_annot.csv', 'condition', 'T2D', 'non-diabetic'),
]:
    expr = pd.read_csv(SC + '\\' + expr_f, index_col=0)
    annot = pd.read_csv(SC + '\\' + annot_f)
    cond = annot[cond_col].values
    # 取可用基因
    avail = [g for g in core_cell_genes if g in expr.index]
    mat = np.vstack([zscore(expr.loc[g].values) for g in avail])
    score = mat.mean(axis=0)
    t2d = cond == t2d_vals
    nd = cond == nd_vals
    p = mw_p(score[t2d], score[nd])
    results[f'{name}_core_score'] = (p, score[nd].mean(), score[t2d].mean(), avail)
    print(f'[{name}] 精简协调评分({"+".join(avail)}): p={p:.3e}, ND mean={score[nd].mean():.3f}, T2D mean={score[t2d].mean():.3f}')
    # 也试试 SLC2A2 + IAPP 双基因
    avail2 = [g for g in ['SLC2A2','IAPP'] if g in expr.index]
    mat2 = np.vstack([zscore(expr.loc[g].values) for g in avail2])
    score2 = mat2.mean(axis=0)
    p2 = mw_p(score2[t2d], score2[nd])
    results[f'{name}_core_score_2'] = (p2, score2[nd].mean(), score2[t2d].mean(), avail2)
    print(f'[{name}] 双基因评分({"+".join(avail2)}): p={p2:.3e}')

# Baron：列名自带 ND_ / T2D_ 前缀
bar = pd.read_csv(SC + r'\baron_beta_expr.csv', index_col=0)
cols = list(bar.columns)
nd_cols = [c for c in cols if c.startswith('ND_')]
t2d_cols = [c for c in cols if c.startswith('T2D_')]
avail = [g for g in core_cell_genes if g in bar.index]
mat = np.vstack([zscore(bar.loc[g].values) for g in avail])
score = mat.mean(axis=0)
t2d_idx = [cols.index(c) for c in t2d_cols]
nd_idx = [cols.index(c) for c in nd_cols]
p = mw_p(score[t2d_idx], score[nd_idx])
results['Baron_core_score'] = (p, score[nd_idx].mean(), score[t2d_idx].mean(), avail)
print(f'[Baron] 精简协调评分({"+".join(avail)}): p={p:.3e}, ND mean={score[nd_idx].mean():.3f}, T2D mean={score[t2d_idx].mean():.3f}')
avail2 = [g for g in ['SLC2A2','IAPP'] if g in bar.index]
mat2 = np.vstack([zscore(bar.loc[g].values) for g in avail2])
score2 = mat2.mean(axis=0)
p2 = mw_p(score2[t2d_idx], score2[nd_idx])
results['Baron_core_score_2'] = (p2, score2[nd_idx].mean(), score2[t2d_idx].mean(), avail2)
print(f'[Baron] 双基因评分({"+".join(avail2)}): p={p2:.3e}')

# ============================================================
# 分析 #2：progenitor/disallowed 评分（bulk）
# 面板 = planA_analysis.R 原始 8 基因（FOXE1,MYCN,CHL1,OPRD1,GRIN2D,ALDH1A3,HES1,SOX9）
# 已验证：R wilcox exact p = 1.369e-11 (GSE76895) / 2.446e-04 (GSE76894)
# ============================================================
prog_genes = ['FOXE1','MYCN','CHL1','OPRD1','GRIN2D','ALDH1A3','HES1','SOX9']

for tag, expr_f, grp_f in [('GSE76895', 'GSE76895_expr.csv', 'GSE76895_groups.csv'),
                            ('GSE76894', 'GSE76894_expr.csv', 'GSE76894_groups.csv')]:
    expr = pd.read_csv(R3 + '\\' + expr_f, index_col=0)
    grp = pd.read_csv(R3 + '\\' + grp_f)
    grp.columns = ['sample', 'group']
    grp = grp[grp['sample'].isin(expr.columns)]
    expr = expr[grp['sample']]
    avail = [g for g in prog_genes if g in expr.index]
    # 用样本内 z-score（跨样本），再跨基因平均
    mat = np.vstack([zscore(expr.loc[g].values) for g in avail])
    score = mat.mean(axis=0)
    t2d = grp['group'].values == 'T2D'
    nd = grp['group'].values == 'ND'
    p = mw_p(score[t2d], score[nd])
    results[f'{tag}_prog_score'] = (p, score[nd].mean(), score[t2d].mean(), avail)
    print(f'[{tag}] progenitor评分({len(avail)} genes: {",".join(avail)}): p={p:.3e}, ND mean={score[nd].mean():.3f}, T2D mean={score[t2d].mean():.3f}')
    # 单基因方向
    for g in avail:
        gp = mw_p(expr.loc[g].values[t2d], expr.loc[g].values[nd])
        print(f'    {g}: T2D mean={expr.loc[g].values[t2d].mean():.3f}, ND mean={expr.loc[g].values[nd].mean():.3f}, p={gp:.3e}')

# ============================================================
# 画图
# ============================================================
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12

# 图 A：单细胞精简协调评分（Xin + Baron 两 panel 合成）
fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))

# Xin core score
xin_expr = pd.read_csv(SC + r'\xin_beta_expr.csv', index_col=0)
xin_annot = pd.read_csv(SC + r'\xin_beta_annot.csv')
xin_cond = xin_annot['condition'].values
xin_avail = [g for g in core_cell_genes if g in xin_expr.index]
xin_mat = np.vstack([zscore(xin_expr.loc[g].values) for g in xin_avail])
xin_score = xin_mat.mean(axis=0)
xin_t2d = xin_cond == 'T2D'
xin_nd = xin_cond == 'non-diabetic'
xin_p = results['Xin_core_score'][0]
ax = axes[0]
bp = ax.boxplot([xin_score[xin_nd], xin_score[xin_t2d]], widths=0.5,
                patch_artist=True, tick_labels=['ND', 'T2DM'])
bp['boxes'][0].set_facecolor('#4C8BF5'); bp['boxes'][1].set_facecolor('#F58A8A')
ax.set_title(f'Xin (GSE81608)\ncore score p={xin_p:.2e}', fontsize=10)
ax.set_ylabel('Core identity score (z)')

# Baron core score
bar_avail = [g for g in core_cell_genes if g in bar.index]
bar_mat = np.vstack([zscore(bar.loc[g].values) for g in bar_avail])
bar_score = bar_mat.mean(axis=0)
bar_p = results['Baron_core_score'][0]
ax = axes[1]
bp = ax.boxplot([bar_score[nd_idx], bar_score[t2d_idx]], widths=0.5,
                patch_artist=True, tick_labels=['ND', 'T2DM'])
bp['boxes'][0].set_facecolor('#4C8BF5'); bp['boxes'][1].set_facecolor('#F58A8A')
ax.set_title(f'Baron (GSE84133)\ncore score p={bar_p:.2e}', fontsize=10)

# progenitor bulk score（发现集 + 验证集 合并成 panel C，双箱）
# 官方 p 值 = R wilcox.test(exact=TRUE) 输出（check_dediff_p.R，与正文一致）
R_P_95 = 1.369e-11
R_P_94 = 2.446e-04
prog_scores = []
prog_ps = [R_P_95, R_P_94]
ax = axes[2]
positions = [1, 2, 3.5, 4.5]
data_groups = []
labels_groups = []
for i, tag in enumerate(['GSE76895', 'GSE76894']):
    expr = pd.read_csv(R3 + '\\' + ('GSE76895_expr.csv' if tag=='GSE76895' else 'GSE76894_expr.csv'), index_col=0)
    grp = pd.read_csv(R3 + '\\' + ('GSE76895_groups.csv' if tag=='GSE76895' else 'GSE76894_groups.csv'))
    grp.columns = ['sample','group']
    grp = grp[grp['sample'].isin(expr.columns)]
    expr = expr[grp['sample']]
    avail = [g for g in prog_genes if g in expr.index]
    mat = np.vstack([zscore(expr.loc[g].values) for g in avail])
    sc = mat.mean(axis=0)
    t2d = grp['group'].values=='T2D'; nd = grp['group'].values=='ND'
    data_groups.append(sc[nd]); data_groups.append(sc[t2d])
    labels_groups.append(f'{tag}\nND'); labels_groups.append(f'{tag}\nT2DM')
bp = ax.boxplot(data_groups, positions=[1,2,3.5,4.5], widths=0.5, patch_artist=True,
                tick_labels=labels_groups)
for j, box in enumerate(bp['boxes']):
    box.set_facecolor('#4C8BF5' if j % 2 == 0 else '#F58A8A')
ax.axvline(2.75, color='gray', linestyle='--', linewidth=0.8)
ax.set_title('Progenitor/disallowed score (bulk)', fontsize=10)
ax.text(0.97, 0.06, f'p = {prog_ps[0]:.2e}\np = {prog_ps[1]:.2e}',
        transform=ax.transAxes, fontsize=8, va='bottom', ha='right', color='#333333')
ax.tick_params(axis='x', labelsize=8, rotation=20)

fig.suptitle('Supplementary: cell-autonomous core score & progenitor-score reversal', y=1.02, fontsize=13)
plt.tight_layout()
fig.savefig(OUT + r'\planA_reviewer_supplement.png', dpi=200, bbox_inches='tight')
fig.savefig(OUT + r'\planA_reviewer_supplement.pdf', bbox_inches='tight')
print('\n图已保存: planA_reviewer_supplement.png/pdf')

# 保存数值结果
# 注意：progenitor/disallowed 评分的 p 值此处为 Python scipy.stats.mannwhitneyu 默认渐近法输出，
# 与稿件正文/图注/Table S8 采用的 R wilcox.test(exact=TRUE) 精确法 p 值存在口径差异：
#   GSE76895：本文件(渐近) 1.44e-09  vs  稿件(R exact) 1.37e-11（check_dediff_p.R 实跑确认）
#   GSE76894：本文件(渐近) 3.73e-04  vs  稿件(R exact) 2.45e-04
# 稿件统一以 R wilcox.test(exact=TRUE) 为最终口径（见 Methods Statistics 段）。
# 核心评分（Xin/Baron）的 p 值两口径一致，无需区分。
with open(BASE + r'\reviewer_supplement_results.txt', 'w', encoding='utf-8') as f:
    f.write('# 口径说明：progenitor/disallowed 评分的 p 值以 R wilcox.test(exact=TRUE) 为最终口径\n')
    f.write('#   GSE76895: R exact p=1.369e-11（正文用）；下方 Python 渐近 p=1.44e-09 仅作参考\n')
    f.write('#   GSE76894: R exact p=2.446e-04（正文用）；下方 Python 渐近 p=3.73e-04 仅作参考\n')
    f.write('# 核心评分（Xin/Baron）Python 与 R 口径一致。\n')
    for k, v in results.items():
        if isinstance(v, tuple):
            f.write(f'{k}: p={v[0]:.6e}, ND_mean={v[1]:.4f}, T2D_mean={v[2]:.4f}, genes={v[3]}\n')
        else:
            f.write(f'{k}: {v}\n')
print('数值已保存: reviewer_supplement_results.txt')
