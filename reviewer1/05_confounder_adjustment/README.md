# Reviewer 1, Comment 5 — Confounder adjustment

## Reviewer comment

> No confounder adjustment. Age, sex, BMI, organ-donor versus pancreatectomy status (IMIDIA contains both), cold ischaemia time and islet purity are all unaddressed. Exocrine contamination alone can move several of these genes.

## Contents

| File | Description |
|---|---|
| `confounder_adjustment.py` | Confounder-adjustment script (ordinary least squares on public data) |
| `confounder_adjustment_results.txt` | Summary of results |
| `confounder_adjustment_output.txt` | Full run output |

## Usage

```bash
export T2D_DATA=/path/to/workspace   # directory containing geo_data/ and r3_inputs/
python confounder_adjustment.py
```

Dependencies: Python standard library + numpy only. The adjustment uses
OLS normal equations (equivalent to
`R lm(score ~ group + age + gender + bmi + batch)`), with p-values from a
Student-t distribution (regularized incomplete beta function),
matching scipy/R.

## Results

> Note: this script uses the same **15-gene uniform cross-cohort panel** (16
> genes minus PDX1, because PDX1 has no probe on the GPL570 platform of the
> validation cohort GSE76894). PDX1 is still analyzed separately as a single
> gene in the discovery cohort.

**The reviewer's concern is confirmed: the identity-score decrease in the
discovery cohort is largely driven by age and is no longer significant after
adjustment.**

| Cohort | Donor source | Samples | identity score unadjusted p | adjusted p | Outcome |
|---|---|---|---|---|---|
| GSE76895 (discovery) | pancreatectomy patients | 68 (36 T2D / 32 ND) | 2.21e-3 | **0.29** | loses significance |
| GSE76894 (validation) | organ donors | 103 (19 T2D / 84 ND) | 3.63e-6 | 6.30e-3 | remains significant |

Key points:

1. **Age is a strong confounder**: in GSE76895 the T2D group averages 65.8 y vs
   60.0 y in ND (a 5.8 y difference); in GSE76894 the T2D group averages 72.2 y
   vs 60.0 y in ND (a 12.2 y difference). The T2D groups are markedly older.
2. **After adjusting for age/gender/bmi/batch**: the identity-score decrease in
   the discovery cohort is no longer significant (p = 0.29); the validation
   cohort remains significant (p = 6.30e-3).
3. Single genes: SLC2A2, PDX1, and NKX6-1 remain significantly down-regulated
   after adjustment (robust); INS, IAPP, and SLC30A8 lose significance in the
   validation cohort after adjustment.
4. **Cold ischaemia time and islet purity**: these two fields are absent from
   the GEO metadata and cannot be adjusted for (stated in the response letter).

## Implications for the manuscript

- **The main conclusion (consistent identity-score decrease) is weakened**: the
  discovery cohort is not significant after adjustment, so the claim of a
  "consistent decrease in both cohorts" should be revised to "robust in the
  validation cohort, directionally consistent but not significant after
  adjustment in the discovery cohort".
- **The single-gene core result (SLC2A2 as the top biomarker) is unaffected**:
  SLC2A2 remains highly significant in both cohorts after adjustment
  (p = 6.0e-9 and 1.1e-5).
- **Overall**: this is a weakening at the level of causal language and strength
  of claims, not a fatal issue, but it should be disclosed explicitly in the
  manuscript.
