# -*- coding: utf-8 -*-
"""
Overfitting assessment of the diagnostic models.

The core concern: the multivariable model's validation AUC of approximately
0.76 is lower than the SLC2A2 single-gene 0.869, with a 19:84 class imbalance,
consistent with overfitting.

This script computes:
1. 5-fold StratifiedKFold cross-validation on the training set (GSE76895, n=68) → generalization AUC
2. Validation-set (GSE76894, n=103) SLC2A2 single-gene vs multivariable logistic AUC
3. DeLong test: whether the multivariable model is significantly better than the SLC2A2 single gene
4. Bootstrap 95% CI (AUC uncertainty)
5. Balanced accuracy + MCC (robust metrics under class imbalance)
6. Training-validation gap (quantification of overfitting)

Data: r3_inputs/GSE76895_expr.csv, GSE76895_groups.csv, GSE76894_expr.csv, GSE76894_groups.csv
Usage: set T2D_DATA to the workspace root, then run.
"""
import os, csv
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, matthews_corrcoef

# Data location via environment variable; by default walks upward to the project root containing r3_inputs/.
def _find_data_root():
    env = os.environ.get('T2D_DATA', '').rstrip('/\\')
    if env and os.path.isdir(env):
        return env
    d = os.path.abspath(os.path.dirname(__file__))
    for _ in range(12):
        if os.path.isdir(os.path.join(d, 'r3_inputs')):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return os.getcwd()

BASE = _find_data_root()
R3IN = os.path.join(BASE, 'r3_inputs')
OUT = os.path.join(BASE, 'ml_output')
os.makedirs(OUT, exist_ok=True)

# 16-gene identity score
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

print('Samples: train(GSE76895)=%d, validation(GSE76894)=%d' % (len(train_y), len(val_y)))
print('Train classes: T2D=%d, ND=%d' % (int(train_y.sum()), int((1-train_y).sum())))
print('Validation classes: T2D=%d, ND=%d' % (int(val_y.sum()), int((1-val_y).sum())))
print()

# common genes (PDX1 has no probe in GSE76894)
common = [g for g in IDG if g in val_expr.index]
print('Common genes: %d (PDX1 has no probe in the validation cohort, excluded)' % len(common))
train_X = train_expr.loc[common].T.values
val_X = val_expr.loc[common].T.values

scaler = StandardScaler().fit(train_X)
Xs_tr = scaler.transform(train_X)
Xs_va = scaler.transform(val_X)

def new_logreg():
    return LogisticRegression(max_iter=2000, C=1.0, random_state=42)

print('=' * 70)
print('1) 5-fold cross-validation on the training set (generalization performance)')
print('=' * 70)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
models = {
    'Logistic regression': new_logreg(),
    'SVM (RBF)': SVC(kernel='rbf', probability=True, C=1.0, gamma='scale', random_state=42),
    'Random forest': RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1),
}
def clone(name):
    if name == 'Logistic regression':
        return new_logreg()
    if name == 'SVM (RBF)':
        return SVC(kernel='rbf', probability=True, C=1.0, gamma='scale', random_state=42)
    return RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)

cv_rows = []
for name in models:
    aucs, baccs, mccs = [], [], []
    for tr_idx, va_idx in cv.split(train_X, train_y):
        s = StandardScaler().fit(train_X[tr_idx])
        Xtr = s.transform(train_X[tr_idx])
        Xva = s.transform(train_X[va_idx])
        m = clone(name).fit(Xtr, train_y[tr_idx])
        prob = m.predict_proba(Xva)[:, 1]
        pred = (prob >= 0.5).astype(int)
        aucs.append(roc_auc_score(train_y[va_idx], prob))
        baccs.append(balanced_accuracy_score(train_y[va_idx], pred))
        mccs.append(matthews_corrcoef(train_y[va_idx], pred))
    print(f'{name}: CV_AUC={np.mean(aucs):.3f}±{np.std(aucs):.3f}, '
          f'balanced_acc={np.mean(baccs):.3f}, MCC={np.mean(mccs):.3f}')
    cv_rows.append([name, f'{np.mean(aucs):.4f}', f'{np.std(aucs):.4f}',
                    f'{np.mean(baccs):.4f}', f'{np.mean(mccs):.4f}'])

print()
print('=' * 70)
print('2) Validation set: SLC2A2 single-gene vs multivariable logistic AUC')
print('=' * 70)
slc_idx = common.index('SLC2A2')
slc_score = -Xs_va[:, slc_idx]  # high expression = ND, low expression = T2D, so negate
slc_auc = roc_auc_score(val_y, slc_score)

logreg = new_logreg().fit(Xs_tr, train_y)
mv_prob = logreg.predict_proba(Xs_va)[:, 1]
mv_auc = roc_auc_score(val_y, mv_prob)

print(f'SLC2A2 single-gene validation AUC = {slc_auc:.3f}')
print(f'Multivariable logistic validation AUC = {mv_auc:.3f}')

# in-sample AUC on the training set (used to quantify the overfitting gap)
slc_train_score = -Xs_tr[:, slc_idx]
slc_train_auc = roc_auc_score(train_y, slc_train_score)
mv_train_prob = logreg.predict_proba(Xs_tr)[:, 1]
mv_train_auc = roc_auc_score(train_y, mv_train_prob)
print()
print(f'SLC2A2 train in-sample AUC = {slc_train_auc:.3f} (train-validation gap {slc_train_auc - slc_auc:+.3f})')
print(f'Multivariable train in-sample AUC = {mv_train_auc:.3f} (train-validation gap {mv_train_auc - mv_auc:+.3f})')

print()
print('=' * 70)
print('3) DeLong test: multivariable vs SLC2A2 single gene')
print('=' * 70)
# DeLong test (standard formula from DeLong et al. 1988)
def delong_test(y, score1, score2):
    """
    DeLong's nonparametric test for the difference between two AUCs.
    Returns (AUC1, AUC2, Z, p). Uses the structural (X, Y) decomposition
    following DeLong et al. (1988), Biometrics 44:837-845.

    V10 is the component attached to each NEGATIVE sample (i = 1..m):
        V10_i = (1/n) * sum_{j in pos} psi(score_pos_j, score_neg_i)
    V01 is the component attached to each POSITIVE sample (j = 1..n):
        V01_j = (1/m) * sum_{i in neg} psi(score_pos_j, score_neg_i)
    where psi(a, b) = 1 if a > b, 0.5 if a == b, 0 otherwise.
    Var(AUC) = S10/m + S01/n.
    """
    from scipy import stats
    score1 = np.asarray(score1, dtype=float)
    score2 = np.asarray(score2, dtype=float)
    y = np.asarray(y)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    m = len(neg)  # number of negatives
    n = len(pos)  # number of positives

    def psi(a, b):
        return (a > b) + 0.5 * (a == b)

    def components(score):
        # V10 (negative samples) and V01 (positive samples)
        V10 = np.empty(m)
        for i, neg_idx in enumerate(neg):
            V10[i] = np.mean(psi(score[pos], score[neg_idx]))
        V01 = np.empty(n)
        for j, pos_idx in enumerate(pos):
            V01[j] = np.mean(psi(score[pos_idx], score[neg]))
        return V10, V01

    V10_1, V01_1 = components(score1)
    V10_2, V01_2 = components(score2)

    # AUC via Mann-Whitney form
    auc1 = np.mean(psi(score1[pos], score1[neg].reshape(-1, 1)))
    auc2 = np.mean(psi(score2[pos], score2[neg].reshape(-1, 1)))

    def var_components(V10a, V01a, V10b, V01b):
        S10 = np.cov(np.vstack([V10a, V10b]))
        S01 = np.cov(np.vstack([V01a, V01b]))
        # Var(AUC_a - AUC_b) = (S10[0,0]+S10[1,1]-2*S10[0,1])/m
        #                        + (S01[0,0]+S01[1,1]-2*S01[0,1])/n
        return (S10[0, 0] + S10[1, 1] - 2 * S10[0, 1]) / m + \
               (S01[0, 0] + S01[1, 1] - 2 * S01[0, 1]) / n

    var_diff = var_components(V10_1, V01_1, V10_2, V01_2)
    z = (auc1 - auc2) / np.sqrt(max(var_diff, 1e-12))
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return auc1, auc2, z, p

auc_mv_dl, auc_slc_dl, z, p_delong = delong_test(val_y, mv_prob, slc_score)
print(f'Multivariable AUC (DeLong) = {auc_mv_dl:.3f}, SLC2A2 AUC (DeLong) = {auc_slc_dl:.3f}')
print(f'DeLong Z = {z:.3f}, p = {p_delong:.4f}')
if p_delong < 0.05 and auc_mv_dl > auc_slc_dl:
    verdict = 'Multivariable significantly better than SLC2A2 single gene'
elif p_delong < 0.05 and auc_mv_dl < auc_slc_dl:
    verdict = 'Multivariable significantly worse than SLC2A2 single gene'
else:
    verdict = 'Multivariable not significantly different from SLC2A2 single gene'
print(f'Conclusion: {verdict}')

print()
print('=' * 70)
print('4) bootstrap 95% CI')
print('=' * 70)
rng = np.random.default_rng(42)
n = len(val_y)
boot_slc, boot_mv = [], []
for _ in range(2000):
    idx = rng.integers(0, n, n)
    if len(np.unique(val_y[idx])) < 2:
        continue
    boot_slc.append(roc_auc_score(val_y[idx], slc_score[idx]))
    boot_mv.append(roc_auc_score(val_y[idx], mv_prob[idx]))
boot_slc, boot_mv = np.array(boot_slc), np.array(boot_mv)
slc_lo, slc_hi = np.percentile(boot_slc, [2.5, 97.5])
mv_lo, mv_hi = np.percentile(boot_mv, [2.5, 97.5])
print(f'SLC2A2 AUC 95% CI = [{slc_lo:.3f}, {slc_hi:.3f}]')
print(f'Multivariable AUC 95% CI = [{mv_lo:.3f}, {mv_hi:.3f}]')

print()
print('=' * 70)
print('5) Class-imbalance-robust metrics (validation set)')
print('=' * 70)
for name in models:
    m = clone(name).fit(Xs_tr, train_y)
    prob = m.predict_proba(Xs_va)[:, 1]
    pred = (prob >= 0.5).astype(int)
    bacc = balanced_accuracy_score(val_y, pred)
    mcc = matthews_corrcoef(val_y, pred)
    print(f'{name}: balanced_acc={bacc:.3f}, MCC={mcc:.3f}')

# save results
with open(os.path.join(OUT, 'r6_overfitting_results.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['metric', 'value'])
    w.writerow(['SLC2A2_val_AUC', f'{slc_auc:.4f}'])
    w.writerow(['SLC2A2_val_CI_lower', f'{slc_lo:.4f}'])
    w.writerow(['SLC2A2_val_CI_upper', f'{slc_hi:.4f}'])
    w.writerow(['multivariable_val_AUC', f'{mv_auc:.4f}'])
    w.writerow(['multivariable_val_CI_lower', f'{mv_lo:.4f}'])
    w.writerow(['multivariable_val_CI_upper', f'{mv_hi:.4f}'])
    w.writerow(['SLC2A2_train_AUC', f'{slc_train_auc:.4f}'])
    w.writerow(['multivariable_train_AUC', f'{mv_train_auc:.4f}'])
    w.writerow(['DeLong_Z', f'{z:.4f}'])
    w.writerow(['DeLong_p', f'{p_delong:.4f}'])

print()
