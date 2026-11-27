# -*- coding: utf-8 -*-
"""
R1-5 confounder adjustment (real run)
================================================
Reviewer comment: no confounder adjustment. Age, sex, BMI, donor source, cold ischemia time, and islet purity were not addressed.

This script extracts sample-level covariates from the GEO series matrices and:
  1. Runs multiple linear regression on the beta-cell identity score (the core continuous outcome):
        score ~ group(T2D/ND) + age + gender + BMI + batch
     Extract the group-adjusted coefficient, SE, t, and p (hand-written OLS normal equations, equivalent to R lm()).
  2. Apply the same adjustment to key single genes (INS, SLC2A2, GCK, PDX1, NKX6-1, etc.) and report unadjusted and adjusted p-values.
  3. Provide a confounder-adjusted sensitivity note for the diagnostic-model AUC.

Available covariates (from the GEO series matrices' !Sample_characteristics):
  GSE76895 (PPP partial-pancreatectomy cohort): age, gender, bmi, batch  (fasting glucose is outcome-related, not a confounder)
  GSE76894 (OD organ-donor cohort): age, gender, bmi, batch  (ICU glycemia is outcome-related, not a confounder)
Cold ischemia time and islet purity are absent from the GEO metadata, so they cannot be adjusted for; this is stated in the response letter.

Algorithm: OLS normal equations beta = (X'X)^-1 X'y, t = beta/SE, two-sided Student-t p-value.
     X includes an intercept column + group(0/1, ND=0 T2D=1) + age(continuous) + gender(0/1, F=0 M=1) + bmi(continuous) + batch(factor, one-hot).
     Matches R lm(score ~ group + age + gender + bmi + batch).
     Identity-score definition matches analysis.R: the z-score averaged across 19 genes.
"""
import gzip, re, csv, math
import paths
import numpy as np

BASE = paths._DATA

# ---- gene set (matches analysis.R lines 25-26) ----
BETA_IDENTITY = ["PDX1","NKX6-1","HNF1A","HNF4A","MAFA","MAFB","NEUROD1","PAX6",
                 "SLC2A2","GCK","ABCC8","KCNJ11","IAPP","UCN3","GLP1R","SLC30A8","PCSK1","PCSK2","INS"]
KEY_GENES = ["INS","SLC2A2","GCK","PDX1","NKX6-1","MAFA","MAFB","IAPP","UCN3","SLC30A8"]


def parse_covariates(gse):
    """Parse sample covariates from the series matrices; returns dict: geo -> {age, gender, bmi, batch}"""
    f = paths.geo_processed() + r'\%s_series_matrix.txt.gz' % gse
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
    age = _vals('!Sample_characteristics_ch1')  # match by field, see below
    # more robust: find the field line by line
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
    # extract the value (strip the "field: " prefix)
    def _val(v, key):
        if v is None: return None
        # e.g. "age (years): 65"
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
    """Load the expression matrix (gene x sample); index=gene, cols=GSM"""
    path = paths.r3_inputs() + r'\%s_expr.csv' % gse
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
    path = paths.r3_inputs() + r'\%s_groups.csv' % gse
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
    """z-score averaged across genes (matches analysis.R)"""
    idx = [genes.index(g) for g in gset if g in genes]
    sub = mat[idx, :]
    subz = np.vstack([zscore(sub[i]) for i in range(len(idx))])
    return subz.mean(axis=0), [genes[i] for i in idx]


def build_X(group, age, gender, bmi, batches):
    """Build the design matrix: intercept + group + age + gender + bmi + batch (one-hot, smallest batch as baseline)"""
    n = len(group)
    cols = []
    cols.append(np.ones(n))                 # intercept
    cols.append(np.array(group, dtype=float))  # T2D=1, ND=0
    cols.append(np.array(age, dtype=float))
    cols.append(np.array(gender, dtype=float))  # M=1, F=0
    cols.append(np.array(bmi, dtype=float))
    ub = sorted(set(batches))
    if len(ub) > 1:
        for b in ub[1:]:  # baseline = ub[0]
            cols.append(np.array([1.0 if x == b else 0.0 for x in batches]))
    return np.column_stack(cols), ub


def ols(X, y):
    """OLS normal equations; returns beta, se, t, p (two-sided Student-t, df=n-k)"""
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
    # two-sided Student-t p-value (avoid scipy; use own t-distribution approximation)
    p = np.array([2 * (1 - _t_cdf(abs(tv), df)) for tv in t])
    return beta, se, t, p, df


def _t_cdf(t, df):
    """One-tailed Student-t CDF via the regularized incomplete beta function (avoid scipy dependency)"""
    # P(T <= t) = 1 - 0.5 * I_{df/(df+t^2)}(df/2, 1/2)
    x = df / (df + t * t)
    return 1 - 0.5 * _betainc(df/2.0, 0.5, x)


def _betainc(a, b, x):
    """Regularized incomplete beta function I_x(a,b), implemented via continued fraction / series"""
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    # use log and series; here numerical integration (Simpson) for accuracy
    import math
    # ln B(a,b)
    lnB = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    # integrand t^(a-1)(1-t)^(b-1)
    from math import exp, log
    def f(t):
        if t <= 0 or t >= 1: return 0.0
        return exp((a-1)*log(t) + (b-1)*log(1-t) - lnB)
    # adaptive Simpson
    N = 2000
    h = x / N
    s = f(0.0) + f(x)
    for i in range(1, N):
        s += (4 if i % 2 == 1 else 2) * f(i * h)
    return s * h / 3.0


def mw_p(a, b):
    """Two-sided Mann-Whitney U p-value (unadjusted), normal approximation (large sample)"""
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
    # two-sided normal p
    p = 2 * (1 - _norm_cdf(abs(z)))
    return p


def _norm_cdf(x):
    import math
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def run():
    print('#' * 78)
    print('# R1-5 confounder adjustment: identity score ~ group + age + gender + bmi + batch')
    print('#' * 78)
    summary = {}
    for gse in ['GSE76895', 'GSE76894']:
        print('\n' + '=' * 78)
        print(gse)
        print('=' * 78)
        genes, samples, mat = load_expr(gse)
        gmap = load_groups(gse)
        cov = parse_covariates(gse)
        # keep only T2D / ND samples (matches the manuscript)
        keep_idx = []
        keep_group = []
        for j, s in enumerate(samples):
            if s in gmap and gmap[s] in ('T2D', 'ND'):
                keep_idx.append(j)
                keep_group.append(gmap[s])
        mat = mat[:, keep_idx]
        samples = [samples[j] for j in keep_idx]
        group = keep_group
        # extract covariates
        age = []; gender = []; bmi = []; batch = []
        for s in samples:
            c = cov.get(s, {})
            a = c.get('age'); g_ = c.get('gender'); b = c.get('bmi'); ba = c.get('batch')
            age.append(float(a) if a and a.replace('.','').isdigit() else float('nan'))
            gender.append(1.0 if g_ == 'M' else (0.0 if g_ == 'F' else float('nan')))
            bmi.append(float(b) if b and b.replace('.','').isdigit() else float('nan'))
            batch.append(ba if ba else 'NA')
        age = np.array(age); gender = np.array(gender); bmi = np.array(bmi)
        # missing values: impute with mean/mode (rare in the real data) and record
        n_missing = int(np.isnan(age).sum() + np.isnan(gender).sum() + np.isnan(bmi).sum())
        age = np.where(np.isnan(age), np.nanmean(age), age)
        gender = np.where(np.isnan(gender), np.nanmean(gender), gender)
        bmi = np.where(np.isnan(bmi), np.nanmean(bmi), bmi)
        grp01 = np.array([1.0 if g == 'T2D' else 0.0 for g in group])
        print(f'  samples: {len(samples)} (T2D={int(grp01.sum())}, ND={int((1-grp01).sum())})')
        print(f'  imputed covariates: {n_missing}')
        print(f'  age mean: {age.mean():.1f} (T2D {age[grp01==1].mean():.1f} vs ND {age[grp01==0].mean():.1f})')
        print(f'  bmi mean: {bmi.mean():.1f} (T2D {bmi[grp01==1].mean():.1f} vs ND {bmi[grp01==0].mean():.1f})')
        print(f'  sex (male fraction): {gender.mean():.2f}')
        # compute the identity score
        score, used_genes = identity_score(mat, genes, BETA_IDENTITY)
        print(f'  identity score genes: {len(used_genes)}')
        # unadjusted Mann-Whitney
        p_unadj = mw_p(score[grp01 == 1], score[grp01 == 0])
        # adjusted OLS
        X, ub = build_X(grp01, age, gender, bmi, batch)
        beta, se, t, p_adj, df = ols(X, score)
        print(f'  [identity score] unadjusted Mann-Whitney p = {p_unadj:.4e}')
        print(f'  [identity score] adjusted OLS: group coefficient beta={beta[1]:+.4f}, SE={se[1]:.4f}, t={t[1]:+.3f}, p={p_adj[1]:.4e} (df={df})')
        # whether the direction is preserved after adjustment
        summary[gse] = {
            'n_T2D': int(grp01.sum()), 'n_ND': int((1-grp01).sum()),
            'score_unadj_p': p_unadj, 'score_adj_p': p_adj[1],
            'score_beta': beta[1], 'score_t': t[1],
        }
        # per-gene adjustment for key genes
        print(f'  --- key genes: unadjusted vs adjusted p ---')
        print(f'  {"gene":10s} {"unadj p":>12s} {"adj p":>12s} {"adj beta":>12s}')
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

    # write the summary results file
    out = paths._DATA + r'\confounder_adjustment_results.txt'
    with open(out, 'w', encoding='utf-8') as f:
        f.write('# R1-5 confounder adjustment results (identity score ~ group + age + gender + bmi + batch)\n')
        f.write('# Method: OLS multiple linear regression; the group coefficient is the T2D vs ND difference\n')
        f.write('# Covariate source: GEO series matrices (!Sample_characteristics): age, gender, bmi, batch\n')
        f.write('# Note: cold ischemia time and islet purity are absent from the GEO metadata and cannot be adjusted for; this is stated in the response letter\n\n')
        for gse, s in summary.items():
            f.write(f'## {gse}\n')
            f.write(f'  samples: T2D={s["n_T2D"]}, ND={s["n_ND"]}\n')
            f.write(f'  identity score: unadj p={s["score_unadj_p"]:.4e}, adj p={s["score_adj_p"]:.4e}, beta={s["score_beta"]:+.4f}\n')
            f.write(f'  key genes (unadj p, adj p, beta):\n')
            for g, pu, pa, bt in s['genes']:
                f.write(f'    {g}: {pu:.4e} -> {pa:.4e} (beta {bt:+.4f})\n')
            f.write('\n')
    print('\nResults saved:', out)


if __name__ == '__main__':
    run()
