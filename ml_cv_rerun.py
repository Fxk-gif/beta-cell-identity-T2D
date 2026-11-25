# -*- coding: utf-8 -*-
"""
多变量诊断模型交叉验证重跑 —— 修复 R1-6/R2-8 空头承诺

背景：回复信承诺"cross-validated performance + AUC 置信区间 + balanced accuracy/MCC"，
但原 ml_identity_model.py 只是简单 train(GSE76895)/test(GSE76894) split，无交叉验证、无 CI。

本脚本：
1. 用 StratifiedKFold 在训练集(GSE76895, n=68)上做 5 折交叉验证，估计三个模型的真实泛化性能
2. 用 bootstrap 估计验证集(GSE76894, n=103)单基因 SLC2A2 AUC 的 95% 置信区间
3. 报告 balanced accuracy 和 Matthews 相关系数(MCC)（类别不平衡 19:84 下的稳健指标）

输出：ml_output/cv_model_metrics.csv（交叉验证指标）、ml_output/auc_ci.csv（AUC 置信区间）
"""
import csv, os
import paths
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, matthews_corrcoef, confusion_matrix

BASE = paths._DATA
R3IN = paths.r3_inputs()
OUT = paths.ml_output()

IDG = ['PDX1','NKX6-1','HNF1A','HNF4A','MAFB','NEUROD1','PAX6',
       'SLC2A2','GCK','ABCC8','IAPP','GLP1R','SLC30A8','PCSK1','PCSK2','INS']

def load(ds):
    df = pd.read_csv(os.path.join(R3IN, f'{ds}_expr.csv'), index_col=0)
    g = pd.read_csv(os.path.join(R3IN, f'{ds}_groups.csv'))
    grp = dict(zip(g['sample'], g['group']))
    samples = [c for c in df.columns]
    y = np.array([1 if grp[s] == 'T2D' else 0 for s in samples])
    return df, y

train_expr, train_y = load('GSE76895')
val_expr, val_y = load('GSE76894')

# 共同基因（PDX1 在 GSE76894 无探针）
common = [g for g in IDG if g in val_expr.index]
train_X = train_expr.loc[common].T.values
val_X = val_expr.loc[common].T.values

# 标准化（训练集拟合）
scaler = StandardScaler().fit(train_X)
Xs_tr = scaler.transform(train_X)
Xs_va = scaler.transform(val_X)

models = {
    'Logistic regression': LogisticRegression(max_iter=2000, C=1.0, random_state=42),
    'SVM (RBF)': SVC(kernel='rbf', probability=True, C=1.0, gamma='scale', random_state=42),
    'Random forest': RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1),
}

def clone_model(name):
    if name == 'Logistic regression':
        return LogisticRegression(max_iter=2000, C=1.0, random_state=42)
    elif name == 'SVM (RBF)':
        return SVC(kernel='rbf', probability=True, C=1.0, gamma='scale', random_state=42)
    else:
        return RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)

print('=' * 70)
print('1) 5 折交叉验证（训练集 GSE76895, n=68, 用 StratifiedKFold）')
print('=' * 70)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_rows = []
for name, mdl in models.items():
    # 交叉验证 AUC（每折在训练折上重新标准化，避免数据泄漏）
    aucs = []
    baccs = []
    mccs = []
    for tr_idx, va_idx in cv.split(train_X, train_y):
        # 折内标准化
        s = StandardScaler().fit(train_X[tr_idx])
        Xtr = s.transform(train_X[tr_idx])
        Xva = s.transform(train_X[va_idx])
        m = clone_model(name)
        m.fit(Xtr, train_y[tr_idx])
        prob = m.predict_proba(Xva)[:, 1]
        pred = (prob >= 0.5).astype(int)
        aucs.append(roc_auc_score(train_y[va_idx], prob))
        baccs.append(balanced_accuracy_score(train_y[va_idx], pred))
        mccs.append(matthews_corrcoef(train_y[va_idx], pred))
    print(f'\n{name}:')
    print(f'  CV AUC (5折, mean±sd): {np.mean(aucs):.3f} ± {np.std(aucs):.3f}')
    print(f'  CV balanced accuracy: {np.mean(baccs):.3f} ± {np.std(baccs):.3f}')
    print(f'  CV MCC: {np.mean(mccs):.3f} ± {np.std(mccs):.3f}')
    cv_rows.append([name, f'{np.mean(aucs):.4f}', f'{np.std(aucs):.4f}',
                    f'{np.mean(baccs):.4f}', f'{np.mean(mccs):.4f}'])

# 保存 CV 结果
with open(os.path.join(OUT, 'cv_model_metrics.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Model', 'CV_AUC_mean', 'CV_AUC_sd', 'CV_balanced_acc_mean', 'CV_MCC_mean'])
    for r in cv_rows:
        w.writerow(r)

print('\n' + '=' * 70)
print('2) 验证集 (GSE76894) 单基因 SLC2A2 AUC bootstrap 置信区间')
print('=' * 70)
# SLC2A2 单基因 AUC + bootstrap CI
slc_idx = common.index('SLC2A2')
slc_va = Xs_va[:, slc_idx]  # 标准化后的 SLC2A2

# 单基因 AUC（用 SLC2A2 表达作为预测得分，注意方向：下调 -> T2D，取负号）
score = -slc_va  # 高表达 = ND，低表达 = T2D，故取负
slc_auc = roc_auc_score(val_y, score)
print(f'SLC2A2 验证集 AUC = {slc_auc:.3f}')

# bootstrap 95% CI（2000 次重采样）
rng = np.random.default_rng(42)
n = len(val_y)
boot_aucs = []
for _ in range(2000):
    idx = rng.integers(0, n, n)
    # 确保两类都有
    if len(np.unique(val_y[idx])) < 2:
        continue
    boot_aucs.append(roc_auc_score(val_y[idx], score[idx]))
boot_aucs = np.array(boot_aucs)
ci_lo, ci_hi = np.percentile(boot_aucs, [2.5, 97.5])
print(f'SLC2A2 AUC 95% CI (bootstrap, 2000次): [{ci_lo:.3f}, {ci_hi:.3f}]')

# 也计算多变量 logistic 的验证集 AUC + CI
logreg = LogisticRegression(max_iter=2000, C=1.0, random_state=42).fit(Xs_tr, train_y)
mv_prob = logreg.predict_proba(Xs_va)[:, 1]
mv_auc = roc_auc_score(val_y, mv_prob)
boot_mv = []
for _ in range(2000):
    idx = rng.integers(0, n, n)
    if len(np.unique(val_y[idx])) < 2:
        continue
    boot_mv.append(roc_auc_score(val_y[idx], mv_prob[idx]))
boot_mv = np.array(boot_mv)
mv_lo, mv_hi = np.percentile(boot_mv, [2.5, 97.5])
print(f'多变量 logistic 验证集 AUC = {mv_auc:.3f}, 95% CI = [{mv_lo:.3f}, {mv_hi:.3f}]')

with open(os.path.join(OUT, 'auc_ci.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['metric', 'AUC', 'CI_lower', 'CI_upper'])
    w.writerow(['SLC2A2_single_gene_val', f'{slc_auc:.4f}', f'{ci_lo:.4f}', f'{ci_hi:.4f}'])
    w.writerow(['multivariable_logistic_val', f'{mv_auc:.4f}', f'{mv_lo:.4f}', f'{mv_hi:.4f}'])

print('\n' + '=' * 70)
print('3) 验证集 balanced accuracy + MCC（类别不平衡稳健指标）')
print('=' * 70)
for name in ['Logistic regression', 'SVM (RBF)', 'Random forest']:
    m = clone_model(name).fit(Xs_tr, train_y)
    prob = m.predict_proba(Xs_va)[:, 1]
    pred = (prob >= 0.5).astype(int)
    bacc = balanced_accuracy_score(val_y, pred)
    mcc = matthews_corrcoef(val_y, pred)
    cm = confusion_matrix(val_y, pred)
    print(f'{name}: balanced_acc={bacc:.3f}, MCC={mcc:.3f}, 混淆矩阵={cm.tolist()}')

print('\n完成。输出: cv_model_metrics.csv, auc_ci.csv')
