# -*- coding: utf-8 -*-
"""
β细胞身份基因诊断预测模型 + 列线图 + 校准曲线
训练集: GSE76895 (68样本: 36 T2D + 32 ND)
验证集: GSE76894 (103样本: 19 T2D + 84 ND, 缺 PDX1)
三模型: Logistic回归 / SVM(RBF) / 随机森林
输出: ROC曲线、混淆矩阵、指标表、列线图、校准曲线
"""
import os, csv, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

warnings.filterwarnings('ignore')

BASE = r'C:/Users/45902/WorkBuddy/2026-09-20-15-12-44'
R3IN = os.path.join(BASE, 'r3_inputs')
OUT  = os.path.join(BASE, 'ml_output')
FIG  = os.path.join(BASE, 'figures', 'png')
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# 中文字体（用于列线图中文标注可选，这里全英文图，避免字体问题）
plt.rcParams['font.family'] = 'DejaVu Sans'

IDG = ['PDX1','NKX6-1','HNF1A','HNF4A','MAFB','NEUROD1','PAX6',
       'SLC2A2','GCK','ABCC8','IAPP','GLP1R','SLC30A8','PCSK1','PCSK2','INS']

def load_expr(ds):
    """读取表达矩阵，返回 (gene行×sample列) DataFrame，行=基因，列=样本"""
    p = os.path.join(R3IN, f'{ds}_expr.csv')
    df = pd.read_csv(p, index_col=0)
    return df  # 行=基因, 列=样本

def load_groups(ds):
    p = os.path.join(R3IN, f'{ds}_groups.csv')
    g = pd.read_csv(p)
    return dict(zip(g['sample'], g['group']))

# ---- 载入训练集 ----
train_expr = load_expr('GSE76895')   # 68 样本
train_grp  = load_groups('GSE76895')
train_samples = [c for c in train_expr.columns]
train_y = np.array([1 if train_grp[s]=='T2D' else 0 for s in train_samples])

# ---- 载入验证集（缺PDX1）----
val_expr = load_expr('GSE76894')
val_grp  = load_groups('GSE76894')
val_samples = [c for c in val_expr.columns]
val_y = np.array([1 if val_grp[s]=='T2D' else 0 for s in val_samples])

# 用两数据集共同基因建模（PDX1 在 GSE76894 无探针，故 15 个共同基因）
common_genes = [g for g in IDG if g in val_expr.index]
train_X = train_expr.loc[common_genes].T.values   # 68 × 15
val_X   = val_expr.loc[common_genes].T.values     # 103 × 15
print(f'训练集: {train_X.shape}, 验证集: {val_X.shape}, 共同基因: {len(common_genes)}')

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score, \
    precision_score, recall_score, f1_score, roc_auc_score

# 标准化（训练集拟合，验证集用同一scaler）
scaler = StandardScaler().fit(train_X)
Xs_tr = scaler.transform(train_X)
Xs_va = scaler.transform(val_X)

models = {
    'Logistic regression': LogisticRegression(max_iter=1000, C=1.0, random_state=42),
    'SVM (RBF)': SVC(kernel='rbf', probability=True, C=1.0, gamma='scale', random_state=42),
    'Random forest': RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1),
}

# ---- 训练 + 预测 ----
results = []
roc_rows = []  # 保存ROC曲线点数据
for name, mdl in models.items():
    mdl.fit(Xs_tr, train_y)
    # 训练集预测概率
    tr_prob = mdl.predict_proba(Xs_tr)[:,1]
    # 验证集预测概率
    va_prob = mdl.predict_proba(Xs_va)[:,1]
    # 验证集 AUC
    va_auc = roc_auc_score(val_y, va_prob)
    # 验证集最佳阈值（约登指数）
    fpr, tpr, th = roc_curve(val_y, va_prob)
    j = tpr - fpr
    best_th = th[np.argmax(j)]
    va_pred = (va_prob >= best_th).astype(int)
    # 指标
    acc = accuracy_score(val_y, va_pred)
    prec = precision_score(val_y, va_pred, zero_division=0)
    rec = recall_score(val_y, va_pred)
    f1 = f1_score(val_y, va_pred, zero_division=0)
    cm = confusion_matrix(val_y, va_pred)
    results.append([name, va_auc, best_th, acc, prec, rec, f1, cm.tolist()])
    roc_rows.append((name, fpr, tpr, va_auc))
    print(f'{name}: 验证集 AUC={va_auc:.3f}, 阈值={best_th:.3f}, ACC={acc:.3f}, '
          f'Precision={prec:.3f}, Recall={rec:.3f}, F1={f1:.3f}')
    print(f'  混淆矩阵(TN,FP/FN,TP): {cm.tolist()}')

# ---- 保存指标表 ----
with open(os.path.join(OUT, 'model_metrics.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['Model', 'Validation_AUC', 'Threshold', 'Accuracy', 'Precision',
                'Recall', 'F1', 'ConfusionMatrix'])
    for r in results:
        w.writerow(r)

# ---- ROC 曲线图 ----
plt.figure(figsize=(7, 6))
colors = {'Logistic regression':'#d62728', 'SVM (RBF)':'#1f77b4', 'Random forest':'#2ca02c'}
for name, fpr, tpr, va_auc in roc_rows:
    plt.plot(fpr, tpr, lw=2.2, color=colors[name],
             label=f'{name} (AUC = {va_auc:.3f})')
plt.plot([0,1],[0,1], 'k--', lw=1.2, label='Chance (AUC = 0.500)')
plt.xlim([-0.02, 1.02]); plt.ylim([-0.02, 1.02])
plt.xlabel('1 - Specificity (False positive rate)', fontsize=12)
plt.ylabel('Sensitivity (True positive rate)', fontsize=12)
plt.title('ROC curves — independent validation set (GSE76894, n = 103)', fontsize=13)
plt.legend(loc='lower right', fontsize=11, frameon=True)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_roc_val.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_roc_val.pdf'))
plt.close()
print('ROC图已保存')

# ---- Logistic 模型系数（供列线图）----
logreg = models['Logistic regression']
coef = logreg.coef_[0]
coef_df = pd.DataFrame({'gene': common_genes, 'coefficient': coef}).sort_values('coefficient', ascending=False)
coef_df.to_csv(os.path.join(OUT, 'logreg_coefficients.csv'), index=False)
print('\nLogistic回归系数(标准化后):')
print(coef_df.to_string(index=False))

# ---- 列线图（Nomogram，手绘基于Logistic系数）----
def draw_nomogram(coef_df, intercept, out_png):
    genes = coef_df['gene'].tolist()
    coefs = coef_df['coefficient'].tolist()
    n = len(genes)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    # 每行基因：标准化表达范围 -2 ~ 2
    xmin, xmax = -2.0, 2.0
    y_top = n + 2.5
    # 分数轴（顶部）
    total_pts = sum(abs(c) for c in coefs) * 2.0
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, y_top)
    for i, (g, c) in enumerate(zip(genes, coefs)):
        y = n - i + 1.5
        # 基因行刻度（标准化表达 -2..2）
        ax.plot([5, 95], [y, y], color='#888888', lw=1.0)
        for v in [-2, -1, 0, 1, 2]:
            x = 5 + (v - xmin)/(xmax-xmin)*90
            ax.plot([x, x], [y-0.12, y+0.12], color='#333333', lw=1.2)
            ax.text(x, y-0.35, f'{v}', ha='center', va='top', fontsize=8)
        ax.text(1.5, y, g, ha='right', va='center', fontsize=9, style='italic')
    # 分数（Points）轴
    ax.plot([5, 95], [n+2.2, n+2.2], color='black', lw=1.5)
    maxpts = max(abs(c) for c in coefs) * 2.0
    for i in range(0, int(maxpts)+1, max(1, int(maxpts//5))):
        x = 5 + i/maxpts*90
        ax.plot([x, x], [n+2.05, n+2.35], color='black', lw=1.0)
        ax.text(x, n+1.85, f'{i}', ha='center', va='top', fontsize=8)
    ax.text(1.5, n+2.2, 'Points', ha='right', va='center', fontsize=10, weight='bold')
    # 总分（Total points）轴
    ax.plot([5, 95], [1.0, 1.0], color='black', lw=1.5)
    tp_range = total_pts
    for i in range(0, int(tp_range)+1, max(1, int(tp_range//10))):
        x = 5 + i/tp_range*90
        ax.plot([x, x], [0.85, 1.15], color='black', lw=1.0)
        ax.text(x, 0.55, f'{i}', ha='center', va='top', fontsize=8)
    ax.text(1.5, 1.0, 'Total Points', ha='right', va='center', fontsize=10, weight='bold')
    # 预测概率轴
    ax.plot([5, 95], [-0.1, -0.1], color='black', lw=1.5)
    probs = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99]
    for p in probs:
        lp = np.log(p/(1-p))
        x = 5 + (lp - intercept) / tp_range * 90
        x = np.clip(x, 5, 95)
        ax.plot([x, x], [-0.25, 0.05], color='black', lw=1.0)
        ax.text(x, -0.55, f'{p:.2f}', ha='center', va='top', fontsize=8)
    ax.text(1.5, -0.1, 'Risk (T2D probability)', ha='right', va='center', fontsize=10, weight='bold')
    ax.axis('off')
    ax.set_title('Nomogram for predicting T2DM based on β-cell identity genes '
                 '(multivariable logistic regression)', fontsize=13)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.savefig(os.path.join(BASE, 'figures', 'ml_nomogram.pdf'))
    plt.close()
    print('列线图已保存')

# 用非标准化系数更直观（这里仍用标准化系数示意，实际投稿可换原始表达）
draw_nomogram(coef_df, logreg.intercept_[0], os.path.join(FIG, 'ml_nomogram.png'))

# ---- 校准曲线（验证集）----
from sklearn.calibration import calibration_curve
plt.figure(figsize=(6.5, 6))
for name in models:
    mdl = models[name]
    va_prob = mdl.predict_proba(Xs_va)[:,1]
    frac_pos, mean_pred = calibration_curve(val_y, va_prob, n_bins=5, strategy='quantile')
    plt.plot(mean_pred, frac_pos, 'o-', lw=2, color=colors[name], label=name)
plt.plot([0,1],[0,1], 'k--', lw=1.2, label='Perfectly calibrated')
plt.xlim([-0.02, 1.02]); plt.ylim([-0.02, 1.02])
plt.xlabel('Mean predicted probability', fontsize=12)
plt.ylabel('Fraction of positives', fontsize=12)
plt.title('Calibration curves (validation set)', fontsize=13)
plt.legend(loc='upper left', fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ml_calibration.png'), dpi=200)
plt.savefig(os.path.join(BASE, 'figures', 'ml_calibration.pdf'))
plt.close()
print('校准曲线已保存')

# ---- 汇总打印 ----
print('\n===== 建模完成 =====')
print('输出目录:', OUT)
print('指标表: model_metrics.csv')
print('Logistic系数: logreg_coefficients.csv')
