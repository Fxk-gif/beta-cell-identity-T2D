# -*- coding: utf-8 -*-
"""
Consistency analysis of the beta-cell identity network's coordinated downregulation
Figure 1: scatter of identity-gene downregulation consistency across the two bulk datasets (log2FC GSE76895 vs GSE76894)
Figure 2: stratified boxplots of the identity score (T2D vs ND, two datasets side by side)
Figure 3: ranking of each identity gene's contribution (loading) to the identity score
Figure 4: scatter of bulk vs single-cell identity-gene downregulation consistency
Figure 5: heatmap of identity-gene expression consistency across datasets (correlation)
"""
import os, csv, warnings
import paths
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

warnings.filterwarnings('ignore')

BASE = paths._DATA
R3IN = paths.r3_inputs()
SCD  = paths.sc_data()
OUT  = paths.ml_output()
FIG  = paths.figures('png')
os.makedirs(OUT, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'

IDG = ['PDX1','NKX6-1','HNF1A','HNF4A','MAFB','NEUROD1','PAX6',
       'SLC2A2','GCK','ABCC8','IAPP','GLP1R','SLC30A8','PCSK1','PCSK2','INS']

def load(ds):
    e = pd.read_csv(os.path.join(R3IN, f'{ds}_expr.csv'), index_col=0)
    g = pd.read_csv(os.path.join(R3IN, f'{ds}_groups.csv'))
    return e, dict(zip(g['sample'], g['group']))

# ============ 1. log2FC consistency between the two datasets ============
def log2fc(ds):
    e, grp = load(ds)
    t2d = [s for s in e.columns if grp[s]=='T2D']
    nd  = [s for s in e.columns if grp[s]=='ND']
    fc = {}
    for g in IDG:
        if g not in e.index: continue
        m_t = e.loc[g, t2d].astype(float).mean()
        m_n = e.loc[g, nd].astype(float).mean()
        fc[g] = np.log2((m_t + 0.01) / (m_n + 0.01))
    return fc

fc95 = log2fc('GSE76895')
fc94 = log2fc('GSE76894')
common = [g for g in IDG if g in fc95 and g in fc94]
x = np.array([fc95[g] for g in common])
y = np.array([fc94[g] for g in common])
r, p = stats.pearsonr(x, y)

# save the FC table
with open(os.path.join(OUT, 'identity_log2fc_both.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gene','log2FC_GSE76895','log2FC_GSE76894'])
    for g in common:
        w.writerow([g, round(fc95[g],4), round(fc94[g],4)])

plt.figure(figsize=(6.5, 6))
plt.axhline(0, color='grey', lw=0.8, ls='--')
plt.axvline(0, color='grey', lw=0.8, ls='--')
plt.scatter(x, y, s=60, c='#d62728', edgecolor='black', zorder=3)
for g, xi, yi in zip(common, x, y):
    plt.annotate(g, (xi, yi), fontsize=8, ha='left', va='bottom')
z = np.polyfit(x, y, 1)
p_ = np.poly1d(z)
xs = np.linspace(min(x)-0.1, max(x)+0.1, 50)
plt.plot(xs, p_(xs), 'b--', lw=1.2, label=f'fit (r = {r:.2f}, p = {p:.2g})')
plt.xlabel('log2FC (GSE76895, discovery)', fontsize=12)
plt.ylabel('log2FC (GSE76894, validation)', fontsize=12)
plt.title('Cross-dataset concordance of identity-gene downregulation', fontsize=12)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_fc_concordance.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_fc_concordance.pdf'))
plt.close()
print(f'FC consistency: r={r:.3f}, p={p:.2g}')

# ============ 2. stratified identity-score boxplots ============
def identity_score(ds):
    e, grp = load(ds)
    genes = [g for g in IDG if g in e.index]
    sub = e.loc[genes].astype(float)
    z = sub.apply(lambda row: (row - row.mean())/row.std(), axis=1)
    score = z.mean(axis=0)
    out = pd.DataFrame({'sample': score.index, 'score': score.values,
                        'group': [grp[s] for s in score.index]})
    return out

sc95 = identity_score('GSE76895')
sc94 = identity_score('GSE76894')

# boxplots
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, (sc, ds, pval) in zip(axes,
                              [(sc95,'GSE76895','8.5e-4'), (sc94,'GSE76894','1.17e-5')]):
    t2d = sc[sc['group']=='T2D']['score']
    nd = sc[sc['group']=='ND']['score']
    u, pv = stats.mannwhitneyu(t2d, nd, alternative='two-sided')
    bp = ax.boxplot([nd, t2d], tick_labels=['ND', 'T2DM'], patch_artist=True,
                    widths=0.55, showmeans=True)
    for patch, c in zip(bp['boxes'], ['#2ca02c', '#d62728']):
        patch.set_facecolor(c); patch.set_alpha(0.55)
    ax.set_title(f'{ds} (Mann-Whitney p = {pv:.2g})', fontsize=12)
    ax.set_ylabel('β-cell identity score (z-score mean)', fontsize=11)
    ax.grid(alpha=0.3, axis='y')
plt.suptitle('β-cell identity score is reduced in T2DM islets', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_identity_box.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_identity_box.pdf'))
plt.close()
print('Identity-score boxplots saved')

# ============ 3. per-gene contribution (loading) ranking ============
# rank by the average log2FC of the two datasets as the "downregulation magnitude"
contrib = pd.DataFrame({'gene': common, 'fc95': x, 'fc94': y})
contrib['mean_fc'] = contrib[['fc95','fc94']].mean(axis=1)
contrib = contrib.sort_values('mean_fc')
plt.figure(figsize=(8, 5.5))
cols = ['#d62728' if v < 0 else '#2ca02c' for v in contrib['mean_fc']]
plt.barh(contrib['gene'], contrib['mean_fc'], color=cols, edgecolor='black', alpha=0.85)
plt.axvline(0, color='grey', lw=0.8)
plt.xlabel('Mean log2FC (T2DM vs ND)', fontsize=12)
plt.ylabel('Identity gene', fontsize=12)
plt.title('Contribution of each identity gene to the coordinated downregulation', fontsize=12)
plt.grid(alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_loading.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_loading.pdf'))
plt.close()
contrib.to_csv(os.path.join(OUT, 'identity_gene_contribution.csv'), index=False)
print('Contribution-ranking figure saved')

# ============ 4. bulk vs single-cell downregulation consistency ============
xin = pd.read_csv(os.path.join(SCD, 'xin_beta_results.csv'))
baron = pd.read_csv(os.path.join(SCD, 'baron_beta_stats.csv'))

# bulk mean log2FC (two datasets)
bulk_fc = pd.DataFrame({'gene': common, 'bulk_fc': (x + y)/2})

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
for ax, sc_df, title in zip(axes,
                            [xin, baron],
                            ['Xin (GSE81608, Smart-seq2)', 'Baron (GSE84133, inDrop)']):
    m = sc_df.merge(bulk_fc, on='gene')
    xx = m['bulk_fc'].values
    yy = m['log2FC'].values
    rr, pp = stats.pearsonr(xx, yy)
    ax.axhline(0, color='grey', lw=0.8, ls='--')
    ax.axvline(0, color='grey', lw=0.8, ls='--')
    ax.scatter(xx, yy, s=55, c='#1f77b4', edgecolor='black', zorder=3)
    for g, xi, yi in zip(m['gene'], xx, yy):
        ax.annotate(g, (xi, yi), fontsize=7, ha='left', va='bottom')
    zz = np.polyfit(xx, yy, 1); pp_ = np.poly1d(zz)
    xs = np.linspace(min(xx)-0.1, max(xx)+0.1, 50)
    ax.plot(xs, pp_(xs), 'r--', lw=1.2, label=f'r = {rr:.2f}, p = {pp:.2g}')
    ax.set_xlabel('Bulk log2FC (mean, two cohorts)', fontsize=11)
    ax.set_ylabel('Single-cell log2FC', fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
plt.suptitle('Concordance between bulk and single-cell downregulation of identity genes',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_bulk_sc_concordance.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_bulk_sc_concordance.pdf'))
plt.close()
print('bulk vs single-cell consistency figure saved')

print('\n===== consistency analysis complete =====')
