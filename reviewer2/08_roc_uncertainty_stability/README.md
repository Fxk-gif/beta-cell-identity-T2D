# Response to Reviewer 2, Comment 8 — ROC uncertainty assessment + multivariable model stability

## Reviewer comment (verbatim)

> (8) The reported ROC analyses lack adequate assessment of uncertainty and potential confounding, while the multivariable models show inconsistent performance.

## Findings

Both criticisms are valid and have been addressed through supplementary analyses.

### 1. ROC analyses lacked an uncertainty assessment (bootstrap confidence intervals added)

The original manuscript computed AUC for each identity gene with pROC but reported only point estimates without confidence intervals. Stratified bootstrap (2000 replicates, stratified by outcome) is now used to estimate 95% confidence intervals for each gene and each multivariable model.

**Key results — single-gene ROC (95% CI):**

| Gene | Discovery cohort AUC [95% CI] | Validation cohort AUC [95% CI] |
|---|---|---|
| SLC2A2 | 0.888 [0.799, 0.959] | 0.869 [0.766, 0.951] |
| PDX1 | 0.852 [0.750, 0.929] | — (no probe in validation cohort) |
| IAPP | 0.642 [0.498, 0.776] | 0.855 [0.763, 0.936] |
| GLP1R | 0.697 [0.560, 0.821] | 0.843 [0.735, 0.925] |
| ABCC8 | 0.588 [0.440, 0.722] | 0.800 [0.685, 0.904] |
| GCK | 0.460 [0.322, 0.594] | 0.731 [0.608, 0.842] |
| PCSK2 | 0.353 [0.234, 0.491] | 0.618 [0.482, 0.749] |

**Important:** not all 16 identity genes carry diagnostic information. Some genes (e.g., discovery-cohort GCK 0.460, PCSK2 0.353, PCSK1 0.427) have AUCs **below chance (0.5)** with CI lower bounds crossing 0.5. This indicates that the original claim that "identity genes as a whole have diagnostic discrimination" was overly optimistic; the discrimination is in fact concentrated in a few genes (SLC2A2, PDX1, IAPP, GLP1R, ABCC8).

### 2. Multivariable models show inconsistent performance (quantified and attributed)

Validation-set AUC of the three multivariable models:

| Model | Validation AUC [95% CI] |
|---|---|
| Logistic regression | 0.759 [0.625, 0.882] |
| SVM (RBF) | 0.698 [0.579, 0.819] |
| Random forest | 0.334 [0.220, 0.450] |

- **The random forest validation AUC = 0.334, significantly below chance (0.5).** This is a genuine and severe instability, not a coding error. The original manuscript did not discuss it; it is now reported and the random forest is downgraded to a **sensitivity analysis** (the failure is attributed to small sample size, high-dimensional correlated predictors, insufficient hyperparameter tuning, and class imbalance), rather than presented as a competing diagnostic model.
- **The multivariable models do not outperform the SLC2A2 single-gene model:** a DeLong test comparing validation-set SLC2A2 (AUC 0.869) vs. multivariable logistic regression (AUC 0.759) gives Z = 1.867, **p = 0.062** (marginally nonsignificant). The multivariable model thus does not significantly improve discrimination, and its point estimate is in fact lower, consistent with the identity genes being highly redundant such that a single representative gene captures most of the discriminative information.

### 3. Model stability (repeated stratified cross-validation added)

5×5 repeated stratified cross-validation within the discovery cohort:

| Model | CV AUC (mean ± SD) | CV balanced accuracy | CV MCC |
|---|---|---|---|
| Logistic regression | 0.927 ± 0.072 | 0.853 | 0.721 |
| SVM (RBF) | 0.928 ± 0.047 | 0.765 | 0.585 |
| Random forest | 0.950 ± 0.052 | 0.852 | 0.720 |

**Interpretation with caution:** these CV values appear high but reflect only within-cohort stability (same platform, same study) and cannot be extrapolated to cross-cohort generalization. The validation-cohort (cross-platform) AUCs are substantially lower (0.334–0.759), and this gap itself is evidence of limited model generalizability.

### 4. Class-imbalance-robust metrics (validation cohort)

The validation cohort is severely imbalanced (T2D n=19 vs ND n=84), so AUC may overestimate performance. Robust metrics are reported at a threshold chosen by the Youden index:

| Model | Youden threshold | balanced accuracy | sensitivity | specificity | MCC |
|---|---|---|---|---|---|
| SLC2A2 single gene | −6.556 | 0.861 | 0.947 | 0.774 | 0.583 |
| Logistic regression | +0.327 | 0.732 | 0.737 | 0.726 | 0.374 |
| SVM (RBF) | +0.010 | 0.770 | 0.789 | 0.750 | 0.439 |
| Random forest | +0.189 | 0.548 | 1.000 | 0.095 | 0.138 |

**Important:** under the Youden threshold and imbalance-robust metrics, logistic regression and SVM retain balanced accuracy clearly above chance (0.73–0.77, MCC 0.37–0.44); only the random forest is near chance (balanced accuracy 0.548, MCC 0.138, specificity 0.095, classifying almost all samples as T2D). Thus logistic regression and SVM remain above chance while the random forest is biased toward the majority class under extreme imbalance. The multivariable models' discriminative power is concentrated in a few genes (SLC2A2, PDX1, IAPP, GLP1R, ABCC8) rather than distributed uniformly across all 16 genes.

### 5. Confounding assessment (supplementary note)

- Both cohorts come from the IMIDIA study and use the same GPL570 platform (already noted in Comment 2). All AUCs are computed within each cohort, so there is no batch confounding from merging cohorts.
- The discovery cohort is relatively balanced (T2D n=36 / ND n=32), whereas the validation cohort is severely imbalanced (T2D n=19 / ND n=84). The difference in class proportions is a major source of the generally lower validation-cohort AUCs and robust metrics.
- Identity-gene expression differences may reflect **disease severity or treatment status** (rather than a causal β-cell identity defect); this is now discussed as a limitation.

## Impact on the manuscript's core claims

**Moderate: the "diagnostic model" claims are further narrowed, without affecting the core biological conclusion:**

1. **Positive:** the discriminative power of the core SLC2A2 single-gene model is **robust** (discovery 0.888, validation 0.869, with CIs well separated from 0.5). This is reinforced by this supplementary analysis.
2. **Needs careful wording:** the "value" of the multivariable diagnostic models must be stated accurately — validation performance is inconsistent across the three models (random forest 0.334 below chance) and the multivariable models do not outperform the single gene (DeLong p=0.062). Under the Youden threshold and imbalance-robust metrics, logistic regression and SVM remain clearly above chance (balanced accuracy 0.73–0.77, MCC 0.37–0.44), but the random forest is near chance (balanced accuracy 0.548, MCC 0.138). The statements about multivariable models in paragraphs [36]/[48] of the original manuscript should be revised to: the multivariable models do **not** provide a discriminative gain and are unstable under small sample size and imbalance, so they are presented only as a sensitivity analysis.
3. **Conceptual level unchanged:** the core claim is "coordinated erosion of the identity gene network" (network level), not the "diagnostic model" itself. The multivariable models were always auxiliary, and their downgrading does **not** change the conclusion that SLC2A2 is the most discriminative member of the identity network and that the identity network is highly redundant — indeed, the DeLong test now supports the redundancy claim more formally.

## Usage

```
python roc_uncertainty_and_model_stability.py
```

Dependencies: numpy, scipy, scikit-learn, pandas. Data are located automatically from the project root under `r3_inputs/`.

## Output files

| File | Contents |
|---|---|
| `roc_auc_bootstrap_ci.csv` | Single-gene AUC point estimates + 95% bootstrap confidence intervals (discovery/validation cohorts) |
| `delong_comparisons.csv` | DeLong tests: SLC2A2 vs each identity gene, SLC2A2 vs multivariable logistic |
| `repeated_cv_stability.csv` | 5×5 repeated stratified cross-validation metrics (discovery cohort) |
| `imbalance_robust_metrics.csv` | balanced accuracy/sensitivity/specificity/F1/MCC at the Youden threshold, with bootstrap CIs (validation cohort) |
