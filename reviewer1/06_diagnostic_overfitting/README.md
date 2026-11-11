# Reviewer 1, Comment #6 — Overfitting assessment code

## Reviewer comment

> The diagnostic underperformance is rationalised rather than explained. A validation AUC of approximately 0.76 for the multivariable model against 0.869 for SLC2A2 alone, with 19 versus 84 class imbalance, is consistent with overfitting. Describing it as "a hallmark of a tightly co-regulated identity network" is post hoc.

## Code contents

`overfitting_assessment.py` computes the following:

1. 5-fold StratifiedKFold cross-validation on the training set (GSE76895, n=68)
2. Validation-set (GSE76894, n=103) AUC for the SLC2A2 single gene vs the multivariable logistic model
3. DeLong test (whether the two AUCs differ significantly)
4. Bootstrap 95% confidence intervals
5. Balanced accuracy + MCC (robust metrics under class imbalance)
6. Training–validation gap (quantification of overfitting)

## Results

| Metric | SLC2A2 single gene | Multivariable logistic |
|---|---|---|
| Validation AUC | **0.869** (CI 0.750–0.951) | **0.759** (CI 0.614–0.889) |
| Training in-sample AUC | 0.888 | 0.981 |
| Training–validation gap | **+0.019** (minimal overfitting) | **+0.222** (severe overfitting) |

**DeLong test**: Z = −1.866, **p = 0.062** — the multivariable model does not
differ significantly from the SLC2A2 single gene (numerically slightly lower),
i.e., the multivariable model adds no discriminative gain.

**5-fold CV AUC** (training set): Logistic 0.952 / SVM 0.928 / Random forest
0.936 — but the Random forest validation MCC = −0.219 (negative, worse than
chance), a typical small-sample overfitting signature.

## Conclusion

The reviewer's judgment is **fully correct**:
1. The multivariable model's training–validation gap of 0.222 is direct
   evidence of overfitting (SLC2A2 is only 0.019)
2. Under the 19:84 class imbalance, the multivariable model does not differ
   significantly from the single gene (DeLong p=0.062) and adds no
   discriminative gain
3. "tightly co-regulated identity network" is indeed a post hoc explanation,
   which has been withdrawn

**Revision to the manuscript**: the conclusion is changed to "**SLC2A2 is the
more parsimonious and better-validated diagnostic feature; the multivariable
model overfits with too few positive samples (n=19) and is not recommended
without external validation on a larger sample**". This further strengthens the
conclusion of SLC2A2 as a single biomarker.

## Usage

```bash
# Environment variable pointing to the workspace root (containing r3_inputs/)
export T2D_DATA="/path/to/workspace"
python overfitting_assessment.py
```

Dependencies: numpy, pandas, scikit-learn, scipy (no other third-party
dependencies).
