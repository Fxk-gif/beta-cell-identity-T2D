# Response to Reviewer 2, Comment 3 — Cross-cohort correlation reanalysis

## What this folder contains

| File | Description |
|---|---|
| `leave_one_out_correlation.py` | Cross-cohort correlation, leave-one-gene-out, and rank-based (Spearman) correlation for the 15 identity genes common to both bulk cohorts. |
| `corrected_log2fc_table.csv` | Corrected per-gene log2FC table (replaces the original Supplementary Table S3). |
| `README.md` | This document. |

## The reviewer's comment

> The reported cross-cohort correlation (Pearson r = 0.831) appears heavily driven by SLC2A2. Recalculation using Supplementary Table S3 shows that excluding SLC2A2 reduces the correlation to approximately 0.08. This substantially weakens the central claim of coordinated network-wide erosion. Leave-one-gene-out and robust correlation analyses are warranted.

## The methodological error (identified and corrected)

Upon re-examination, we identified an error in the original calculation underlying Supplementary Table S3 and the reported cross-cohort correlation.

The expression matrices (GSE76895 and GSE76894, both GPL570) are RMA-normalized and therefore **already in log2 space**. The original script incorrectly applied an additional `log2(ratio)` transformation to these log2 values (i.e., a **double-log transformation**):

```python
# Original (incorrect):
fc[g] = np.log2((mean_T2D + 0.01) / (mean_ND + 0.01))
```

The correct fold change for log2-space data is the **difference** of group means in log2 space:

```python
# Corrected:
fc[g] = mean_T2D - mean_ND
```

This error compressed the magnitude of the fold changes and, for several genes, reversed their sign.

## Corrected results

### Per-gene corrected log2FC

| Gene | GSE76895 | GSE76894 |
|---|---|---|
| NKX6-1 | −0.4517 | −0.2804 |
| HNF1A | −0.3204 | −0.2017 |
| HNF4A | −0.1426 | +0.1355 |
| MAFB | −0.0653 | −0.2618 |
| NEUROD1 | +0.0580 | −0.5167 |
| PAX6 | +0.0239 | −0.4846 |
| SLC2A2 | −1.4181 | −2.0017 |
| GCK | +0.0624 | −0.5920 |
| ABCC8 | −0.0608 | −0.9815 |
| IAPP | −0.0656 | −0.7878 |
| GLP1R | −0.4091 | −0.7493 |
| SLC30A8 | +0.0049 | −0.4836 |
| PCSK1 | +0.1144 | −0.7142 |
| PCSK2 | +0.1898 | −0.4362 |
| INS | −0.0051 | −0.2538 |

### Cross-cohort correlation (corrected, 15 genes)

| Metric | Value | p-value |
|---|---|---|
| Pearson r | **0.668** | 6.5 × 10⁻³ |
| Spearman ρ | **0.004** | 0.99 |

### Leave-one-gene-out (Pearson)

| Excluded gene | r |
|---|---|
| SLC2A2 | **−0.195** |
| any other single gene | +0.66 to +0.73 |

### Leave-one-gene-out (Spearman, rank-based)

| Excluded gene | ρ |
|---|---|
| SLC2A2 | **−0.226** |
| any other single gene | −0.15 to +0.10 |

### Direction-consistency (14 non-SLC2A2 genes)

Only **7 of 14** genes change in the same direction across the two cohorts; the other 7 change in opposite directions (HNF4A, NEUROD1, PAX6, GCK, SLC30A8, PCSK1, PCSK2).

## Conclusion

Under the corrected calculation, the positive cross-cohort correlation is driven almost entirely by SLC2A2 (the only gene robustly and consistently downregulated in both cohorts, log2FC = −1.42 and −2.00). There is no evidence of coordinated, network-wide directional change across the broader HNF1A–PDX1 / maturity gene set in the bulk data. The corresponding manuscript statements (Abstract, Results, Discussion) and Supplementary Table S3 have been revised accordingly.
