# -*- coding: utf-8 -*-
"""
Supplementary analyses for the reviewer response (revised after the donor-level downgrade).

Background: the initial submission built a "core score" from SLC2A2/IAPP/INS at cell level and
interpreted the significant downregulation as "cell-autonomous" identity loss. After Reviewers R1-1/R2-1 flagged the cell-level
test as pseudoreplication, this repository now re-runs the analysis at donor level (see xin_donor_level.py / baron_donor_level.py);
the conclusion is that no gene is significantly downregulated at donor level (SLC2A2 FDR=0.109, IAPP p=0.64, core score p=0.85),
and the only significant changes are GCK/PCSK2 upregulation (FDR=0.042, opposite to the bulk direction).

This script retains two supplementary analyses (matching the manuscript: bulk uses Wilcoxon / scipy mannwhitneyu, z-score averaged across genes):
#1 single-cell core score (SLC2A2 + IAPP + INS) -- exploratory reference only, noted as not significant at donor level;
#2 progenitor/disallowed gene score (bulk) -- verifies it is lower in T2DM (supporting reduced mature beta-cell program rather than dedifferentiation).
Note: this script is still a cell-level implementation (legacy); its results are NOT used as the basis for any manuscript conclusion;
    the manuscript single-cell conclusions follow the donor-level scripts (xin_donor_level.py / baron_donor_level.py).
"""
import numpy as np
import paths
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = paths._DATA
SC = paths.sc_data()
R3 = paths.r3_inputs()
OUT = paths.figures()
import os
os.makedirs(OUT, exist_ok=True)

def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)  # matches R sd() (n-1)
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)

def mw_p(a, b):
    return stats.mannwhitneyu(a, b, alternative='two-sided').pvalue

results = {}

# ============================================================
# Analysis #1: single-cell core score (SLC2A2 + IAPP + INS) -- exploratory, not significant at donor level
# ============================================================
core_cell_genes = ['SLC2A2', 'IAPP', 'INS']  # genes of the initial "core score" (cell-level, legacy)

for name, expr_f, annot_f, cond_col, t2d_vals, nd_vals in [
    ('Xin', 'xin_beta_expr.csv', 'xin_beta_annot.csv', 'condition', 'T2D', 'non-diabetic'),
]:
    expr = pd.read_csv(SC + '\\' + expr_f, index_col=0)
    annot = pd.read_csv(SC + '\\' + annot_f)
    cond = annot[cond_col].values
    # keep available genes
    avail = [g for g in core_cell_genes if g in expr.index]
    mat = np.vstack([zscore(expr.loc[g].values) for g in avail])
    score = mat.mean(axis=0)
    t2d = cond == t2d_vals
    nd = cond == nd_vals
    p = mw_p(score[t2d], score[nd])
    results[f'{name}_core_score'] = (p, score[nd].mean(), score[t2d].mean(), avail)
    print(f'[{name}] core score({"+".join(avail)}): p={p:.3e}, ND mean={score[nd].mean():.3f}, T2D mean={score[t2d].mean():.3f}')
    # also try the SLC2A2 + IAPP two-gene score
    avail2 = [g for g in ['SLC2A2','IAPP'] if g in expr.index]
    mat2 = np.vstack([zscore(expr.loc[g].values) for g in avail2])
    score2 = mat2.mean(axis=0)
    p2 = mw_p(score2[t2d], score2[nd])
    results[f'{name}_core_score_2'] = (p2, score2[nd].mean(), score2[t2d].mean(), avail2)
    print(f'[{name}] two-gene score({"+".join(avail2)}): p={p2:.3e}')

# Baron: columns carry ND_ / T2D_ prefixes
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
print(f'[Baron] core score({"+".join(avail)}): p={p:.3e}, ND mean={score[nd_idx].mean():.3f}, T2D mean={score[t2d_idx].mean():.3f}')
avail2 = [g for g in ['SLC2A2','IAPP'] if g in bar.index]
mat2 = np.vstack([zscore(bar.loc[g].values) for g in avail2])
score2 = mat2.mean(axis=0)
p2 = mw_p(score2[t2d_idx], score2[nd_idx])
results['Baron_core_score_2'] = (p2, score2[nd_idx].mean(), score2[t2d_idx].mean(), avail2)
print(f'[Baron] two-gene score({"+".join(avail2)}): p={p2:.3e}')

# ============================================================
# Analysis #2: progenitor/disallowed score (bulk)
# panel = the original 8 genes from planA_analysis.R (FOXE1,MYCN,CHL1,OPRD1,GRIN2D,ALDH1A3,HES1,SOX9)
# verified: R wilcox exact p = 1.369e-11 (GSE76895) / 2.446e-04 (GSE76894)
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
    # within-sample z-score (across samples), then averaged across genes
    mat = np.vstack([zscore(expr.loc[g].values) for g in avail])
    score = mat.mean(axis=0)
    t2d = grp['group'].values == 'T2D'
    nd = grp['group'].values == 'ND'
    p = mw_p(score[t2d], score[nd])
    results[f'{tag}_prog_score'] = (p, score[nd].mean(), score[t2d].mean(), avail)
    print(f'[{tag}] progenitor score({len(avail)} genes: {",".join(avail)}): p={p:.3e}, ND mean={score[nd].mean():.3f}, T2D mean={score[t2d].mean():.3f}')
    # per-gene direction
    for g in avail:
        gp = mw_p(expr.loc[g].values[t2d], expr.loc[g].values[nd])
        print(f'    {g}: T2D mean={expr.loc[g].values[t2d].mean():.3f}, ND mean={expr.loc[g].values[nd].mean():.3f}, p={gp:.3e}')

# ============================================================
# plot
# ============================================================
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12

# Panel A: single-cell core score (Xin + Baron combined)
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

# progenitor bulk score (discovery + validation combined into panel C, two boxes)
# official p-values = R wilcox.test(exact=TRUE) output (check_dediff_p.R, matches the manuscript)
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
print('\nFigures saved: planA_reviewer_supplement.png/pdf')

# save the numerical results
# Note: the progenitor/disallowed p-values here are from Python scipy.stats.mannwhitneyu (asymptotic by default),
# which differ from the exact R wilcox.test(exact=TRUE) p-values used in the manuscript text / legend / Table S8:
#   GSE76895: this file (asymptotic) 1.44e-09  vs  manuscript (R exact) 1.37e-11 (confirmed by check_dediff_p.R)
#   GSE76894: this file (asymptotic) 3.73e-04  vs  manuscript (R exact) 2.45e-04
# The manuscript uses R wilcox.test(exact=TRUE) as the final standard (see the Methods, Statistics section).
# The core-score (Xin/Baron) p-values agree between the two methods, so no distinction is needed.
with open(BASE + r'\reviewer_supplement_results.txt', 'w', encoding='utf-8') as f:
    f.write('# Standard: the progenitor/disallowed p-values follow R wilcox.test(exact=TRUE)\n')
    f.write('#   GSE76895: R exact p=1.369e-11 (used in the manuscript); the Python asymptotic p=1.44e-09 below is for reference only\n')
    f.write('#   GSE76894: R exact p=2.446e-04 (used in the manuscript); the Python asymptotic p=3.73e-04 below is for reference only\n')
    f.write('# The core-score (Xin/Baron) p-values agree between Python and R.\n')
    for k, v in results.items():
        if isinstance(v, tuple):
            f.write(f'{k}: p={v[0]:.6e}, ND_mean={v[1]:.4f}, T2D_mean={v[2]:.4f}, genes={v[3]}\n')
        else:
            f.write(f'{k}: {v}\n')
print('Numerical results saved: reviewer_supplement_results.txt')
