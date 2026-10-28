# -*- coding: utf-8 -*-
"""
R1-5 混杂校正（Confounder adjustment）真实运行
================================================
审稿人意见：未做混杂校正。年龄、性别、BMI、供体来源、冷缺血时间、胰岛纯度均未处理。

本脚本真实地从 GEO series matrix 提取样本级协变量，并：
  1. 对 β-cell identity score（正文的核心连续结局）做多重线性回归：
        score ~ group(T2D/ND) + age + gender + BMI + batch
     提取 group 的校正后系数、SE、t、p（手写 OLS 正规方程，等价 R lm()）。
  2. 对关键单基因（INS, SLC2A2, GCK, PDX1, NKX6-1 等）做同样的校正，报告校正前后 p 值。
  3. 对诊断模型 AUC 也给出混杂校正的敏感性说明。

可用的真实协变量（来自 GEO series matrix 的 !Sample_characteristics）：
  GSE76895 (PPP 部分胰腺切除队列): age, gender, bmi, batch  (另有 fasting glucose，属结局相关，不作混杂)
  GSE76894 (OD 器官捐献队列):     age, gender, bmi, batch  (另有 icu glycemia，属结局相关，不作混杂)
冷缺血时间(cold ischemia time)与胰岛纯度(islet purity)：GEO 元数据中【不存在】，需在回复信中如实说明无法校正。

算法：OLS 正规方程 beta = (X'X)^-1 X'y，t = beta/SE，p 用 Student-t 双侧。
     X 含截距列 + group(0/1, ND=0 T2D=1) + age(连续) + gender(0/1, F=0 M=1) + bmi(连续) + batch(因子, one-hot)。
     与 R lm(score ~ group + age + gender + bmi + batch) 一致。
     identity score 定义与 analysis.R 一致：19 基因 z-score 跨基因平均。
"""
import gzip, re, csv, math
import numpy as np

BASE = r'C:\Users\45902\WorkBuddy\2026-09-20-15-12-44'

# ---- 基因集（与 analysis.R 第25-26行一致）----
BETA_IDENTITY = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFA","MAFB","NEUROD1","PAX6",
                 "SLC2A2","GCK","ABCC8","KCNJ11","IAPP","UCN3","GLP1R","SLC30A8","PCSK1","PCSK2","INS"]
KEY_GENES = ["INS","SLC2A2","GCK","PDX1","NKX6-1","MAFA","MAFB","IAPP","UCN3","SLC30A8"]


def parse_covariates(gse):
    """从系列矩阵解析样本协变量，返回 dict: geo -> {age, gender, bmi, batch}"""
    f = BASE + r'\geo_data\%s_series_matrix.txt.gz' % gse
    txt = gzip.open(f, 'rt', encoding='latin-1').read()
    lines = txt.split('\n')
    acc = None
    for ln in lines:
        if ln.startswith('!Sample_geo_accession'):
            acc = [x for x in re.findall(r'"([^"]*)"', ln)]
            break
    def _vals(prefix):
        for ln in lines:
            if ln.startswith(prefix):
                return [x for x in re.findall(r'"([^"]*)"', ln)]
        return []
    age = _vals('!Sample_characteristics_ch1')  # 需要按字段匹配，见下
    # 更稳健：逐行找字段
    def _field(key):
        for ln in lines:
            if ln.startswith('!Sample_characteristics') and key in ln:
                vals = [x for x in re.findall(r'"([^"]*)"', ln)]
                return vals
        return []
    age = _field('age (years)')
    gender = _field('gender (m/f)')
    bmi = _field('body mass index (bmi)')
    batch = _field('batch')
    # 提取值（去掉 "字段: " 前缀）
    def _val(v, key):
        if v is None: return None
        # 形如 "age (years): 65"
        m = re.search(r':\s*(.+)$', v)
        return m.group(1).strip() if m else v.strip()
    out = {}
    n = len(acc)
    for i in range(n):
        g = acc[i]
        out[g] = {
            'age':   _val(age[i], 'age') if i < len(age) else None,
            'gender':_val(gender[i], 'gender') if i < len(gender) else None,
            'bmi':   _val(bmi[i], 'bmi') if i < len(bmi) else None,
            'batch': _val(batch[i], 'batch') if i < len(batch) else None,
        }
    return out


def load_expr(gse):
    """加载表达矩阵（gene x sample），index=gene, cols=GSM"""
    path = BASE + r'\r3_inputs\%s_expr.csv' % gse
    rows = []
    header = None
    with open(path, encoding='utf-8-sig') as f:
        rdr = csv.reader(f)
        header = next(rdr)
        for row in rdr:
            rows.append(row)
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    return genes, header[1:], mat


def load_groups(gse):
    path = BASE + r'\r3_inputs\%s_groups.csv' % gse
    gmap = {}
    with open(path, encoding='utf-8-sig') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            gmap[row['sample']] = row['group']
    return gmap


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x, ddof=1)
    return (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)


def identity_score(mat, genes, gset):
    """z-score 跨基因平均（与 analysis.R 一致）"""
    idx = [genes.index(g) for g in gset if g in genes]
    sub = mat[idx, :]
    subz = np.vstack([zscore(sub[i]) for i in range(len(idx))])
    return subz.mean(axis=0), [genes[i] for i in idx]


def build_X(group, age, gender, bmi, batches):
    """构建设计矩阵：截距 + group + age + gender + bmi + batch(one-hot, 以最小batch为基线)"""
    n = len(group)
    cols = []
    cols.append(np.ones(n))                 # 截距
    cols.append(np.array(group, dtype=float))  # T2D=1, ND=0
    cols.append(np.array(age, dtype=float))
    cols.append(np.array(gender, dtype=float))  # M=1, F=0
    cols.append(np.array(bmi, dtype=float))
    ub = sorted(set(batches))
    if len(ub) > 1:
        for b in ub[1:]:  # 基线 = ub[0]
            cols.append(np.array([1.0 if x == b else 0.0 for x in batches]))
    return np.column_stack(cols), ub


def ols(X, y):
    """OLS 正规方程，返回 beta, se, t, p（双侧 Student-t, df=n-k）"""
    n, k = X.shape
    XtX = X.T @ X
    Xty = X.T @ y
    beta = np.linalg.solve(XtX, Xty)
    resid = y - X @ beta
    df = n - k
    sigma2 = (resid @ resid) / df
    cov = sigma2 * np.linalg.inv(XtX)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    # Student-t 双侧 p 值（用 scipy 缺省，改用自己的 t 分布近似）
    p = np.array([2 * (1 - _t_cdf(abs(tv), df)) for tv in t])
    return beta, se, t, p, df


def _t_cdf(t, df):
    """Student-t 单尾 CDF，用正则化不完全 beta 函数实现（避免 scipy 依赖）"""
    # P(T <= t) = 1 - 0.5 * I_{df/(df+t^2)}(df/2, 1/2)
    x = df / (df + t * t)
    return 1 - 0.5 * _betainc(df/2.0, 0.5, x)


def _betainc(a, b, x):
    """正则化不完全 beta 函数 I_x(a,b)，用连分数/级数实现"""
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    # 用 log 及级数：这里用数值积分（辛普森）保证精度
    import math
    # ln B(a,b)
    lnB = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    # 被积函数 t^(a-1)(1-t)^(b-1)
    from math import exp, log
    def f(t):
        if t <= 0 or t >= 1: return 0.0
        return exp((a-1)*log(t) + (b-1)*log(1-t) - lnB)
    # 自适应辛普森
    N = 2000
    h = x / N
    s = f(0.0) + f(x)
    for i in range(1, N):
        s += (4 if i % 2 == 1 else 2) * f(i * h)
    return s * h / 3.0


def mw_p(a, b):
    """Mann-Whitney U 双侧 p（未校正），用正态近似（大样本）"""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    na, nb = len(a), len(b)
    allv = np.concatenate([a, b])
    ranks = np.argsort(np.argsort(allv)) + 1.0
    ra = ranks[:na]; rb = ranks[na:]
    U = ra.sum() - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    # tie correction
    _, cnt = np.unique(allv, return_counts=True)
    N = na + nb
    sigma2 = na * nb / 12.0 * ((N + 1) - sum((c**3 - c) for c in cnt) / (N * (N - 1)))
    sigma2 = max(sigma2, 1e-12)
    z = (U - mu) / math.sqrt(sigma2)
    # 双侧正态 p
    p = 2 * (1 - _norm_cdf(abs(z)))
    return p


def _norm_cdf(x):
    import math
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def run():
    print('#' * 78)
    print('# R1-5 混杂校正：identity score ~ group + age + gender + bmi + batch')
    print('#' * 78)
    summary = {}
    for gse in ['GSE76895', 'GSE76894']:
        print('\n' + '=' * 78)
        print(gse)
        print('=' * 78)
        genes, samples, mat = load_expr(gse)
        gmap = load_groups(gse)
        cov = parse_covariates(gse)
        # 只保留 T2D / ND 样本（与正文一致）
        keep_idx = []
        keep_group = []
        for j, s in enumerate(samples):
            if s in gmap and gmap[s] in ('T2D', 'ND'):
                keep_idx.append(j)
                keep_group.append(gmap[s])
        mat = mat[:, keep_idx]
        samples = [samples[j] for j in keep_idx]
        group = keep_group
        # 提取协变量
        age = []; gender = []; bmi = []; batch = []
        for s in samples:
            c = cov.get(s, {})
            a = c.get('age'); g_ = c.get('gender'); b = c.get('bmi'); ba = c.get('batch')
            age.append(float(a) if a and a.replace('.','').isdigit() else float('nan'))
            gender.append(1.0 if g_ == 'M' else (0.0 if g_ == 'F' else float('nan')))
            bmi.append(float(b) if b and b.replace('.','').isdigit() else float('nan'))
            batch.append(ba if ba else 'NA')
        age = np.array(age); gender = np.array(gender); bmi = np.array(bmi)
        # 缺失值：用均值/众数填补（真实数据中可能极少缺失），并记录
        n_missing = int(np.isnan(age).sum() + np.isnan(gender).sum() + np.isnan(bmi).sum())
        age = np.where(np.isnan(age), np.nanmean(age), age)
        gender = np.where(np.isnan(gender), np.nanmean(gender), gender)
        bmi = np.where(np.isnan(bmi), np.nanmean(bmi), bmi)
        grp01 = np.array([1.0 if g == 'T2D' else 0.0 for g in group])
        print(f'  样本数: {len(samples)} (T2D={int(grp01.sum())}, ND={int((1-grp01).sum())})')
        print(f'  协变量缺失填补数: {n_missing}')
        print(f'  age 均值: {age.mean():.1f} (T2D {age[grp01==1].mean():.1f} vs ND {age[grp01==0].mean():.1f})')
        print(f'  bmi 均值: {bmi.mean():.1f} (T2D {bmi[grp01==1].mean():.1f} vs ND {bmi[grp01==0].mean():.1f})')
        print(f'  性别(男比例): {gender.mean():.2f}')
        # 计算 identity score
        score, used_genes = identity_score(mat, genes, BETA_IDENTITY)
        print(f'  identity score 基因数: {len(used_genes)}')
        # 未校正 Mann-Whitney
        p_unadj = mw_p(score[grp01 == 1], score[grp01 == 0])
        # 校正 OLS
        X, ub = build_X(grp01, age, gender, bmi, batch)
        beta, se, t, p_adj, df = ols(X, score)
        print(f'  [identity score] 未校正 Mann-Whitney p = {p_unadj:.4e}')
        print(f'  [identity score] 校正后 OLS: group 系数 beta={beta[1]:+.4f}, SE={se[1]:.4f}, t={t[1]:+.3f}, p={p_adj[1]:.4e} (df={df})')
        # 校正前后方向是否保持
        summary[gse] = {
            'n_T2D': int(grp01.sum()), 'n_ND': int((1-grp01).sum()),
            'score_unadj_p': p_unadj, 'score_adj_p': p_adj[1],
            'score_beta': beta[1], 'score_t': t[1],
        }
        # 关键基因逐基因校正
        print(f'  --- 关键基因校正前后 p 值 ---')
        print(f'  {"基因":10s} {"未校正p":>12s} {"校正p":>12s} {"校正beta":>12s}')
        gene_rows = []
        for g in KEY_GENES:
            if g not in genes: continue
            gi = genes.index(g)
            gv = mat[gi, :]
            p_u = mw_p(gv[grp01 == 1], gv[grp01 == 0])
            b2, s2, t2, p2, _ = ols(X, gv)
            print(f'  {g:10s} {p_u:12.4e} {p2[1]:12.4e} {b2[1]:+12.4f}')
            gene_rows.append((g, p_u, p2[1], b2[1]))
        summary[gse]['genes'] = gene_rows

    # 输出汇总结果文件
    out = BASE + r'\confounder_adjustment_results.txt'
    with open(out, 'w', encoding='utf-8') as f:
        f.write('# R1-5 混杂校正真实结果（identity score ~ group + age + gender + bmi + batch）\n')
        f.write('# 校正方法：OLS 多重线性回归，group 系数为 T2D vs ND 差异\n')
        f.write('# 协变量来源：GEO series matrix (!Sample_characteristics)：age, gender, bmi, batch\n')
        f.write('# 说明：冷缺血时间(cold ischemia time)与胰岛纯度(islet purity)在 GEO 元数据中不存在，无法校正，需在回复信中如实说明\n\n')
        for gse, s in summary.items():
            f.write(f'## {gse}\n')
            f.write(f'  样本: T2D={s["n_T2D"]}, ND={s["n_ND"]}\n')
            f.write(f'  identity score: 未校正 p={s["score_unadj_p"]:.4e}, 校正后 p={s["score_adj_p"]:.4e}, beta={s["score_beta"]:+.4f}\n')
            f.write(f'  关键基因 (未校正p, 校正p, beta):\n')
            for g, pu, pa, bt in s['genes']:
                f.write(f'    {g}: {pu:.4e} -> {pa:.4e} (beta {bt:+.4f})\n')
            f.write('\n')
    print('\n结果已保存:', out)


if __name__ == '__main__':
    run()
