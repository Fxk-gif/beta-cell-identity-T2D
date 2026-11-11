# Response to Reviewer 2, Comment 5 — Probe annotation and uniform gene panel

## What this folder contains

| File | Description |
|---|---|
| `uniform_panel_reanalysis.py` | Probe-annotation check and reanalysis of the identity score using a uniform 15-gene panel across both cohorts. |
| `uniform_15gene_identity_scores.csv` | Per-sample identity score computed with the uniform 15-gene panel. |
| `README.md` | This document. |

## The reviewer's comment

> The manuscript describes a 16-gene identity score, although PDX1 is reportedly unavailable in the validation cohort. Moreover, PDX1 coverage is inconsistent between Figure 3 and Supplementary Table S5. The authors must clarify probe annotations, reconcile these discrepancies, and repeat the analyses using an identical gene panel across cohorts.

## 1. PDX1 probe annotation (clarified)

The GPL570 platform (Affymetrix HG-U133 Plus 2.0) contains three PDX1 probes: `208559_at`, `210937_s_at`, and `210938_at`. Their presence in the raw series-matrix files is:

| Probe | GSE76895 (discovery) | GSE76894 (validation) |
|---|---|---|
| 208559_at | absent | absent |
| 210937_s_at | **present** | absent |
| 210938_at | **present** | absent |

Therefore PDX1 is genuinely unmeasurable in the validation cohort (GSE76894): all three PDX1 probes are absent from its raw matrix. This is a property of the data, not an annotation error.

## 2. Figure 3 vs Supplementary Table S5 (reconciled)

The inconsistency arose because:
- **Figure 3** (single-gene ROC diagnostic performance) presented PDX1, but this was based on the **discovery cohort only**.
- **Supplementary Table S5** (regulatory-network node table) left the PDX1 log2FC columns blank ("—") because PDX1 could not be quantified in the validation cohort.

The two are now made consistent: Figure 3 explicitly states the PDX1 ROC is discovery-cohort only, and Table S5 carries a footnote explaining that PDX1 log2FC is not reported because PDX1 was not measurable in the validation cohort.

## 3. Uniform 15-gene panel reanalysis

Because PDX1 is unavailable in the validation cohort, the only way to obtain an identical panel across cohorts is to remove PDX1. We recomputed the β-cell identity score using the uniform **15-gene panel** (16 genes minus PDX1) in both cohorts, using the identical z-score-averaging formula.

### Results

| Cohort | T2D mean | ND mean | Mann-Whitney U p |
|---|---|---|---|
| GSE76895 (discovery, n=68) | −0.1023 | +0.1151 | **2.21 × 10⁻³** |
| GSE76894 (validation, n=103) | −0.6911 | +0.1563 | **3.63 × 10⁻⁶** |

### Impact of dropping PDX1 (discovery cohort)

The 16-gene score and the 15-gene score are nearly identical in the discovery cohort (Pearson r = 0.9965; 16-gene p = 8.5 × 10⁻⁴ vs 15-gene p = 2.21 × 10⁻³). Removing PDX1 therefore does not materially change the identity-score conclusion.

## Conclusion

Using a fully identical 15-gene panel across both cohorts, the β-cell identity score is significantly reduced in T2D in both the discovery (p = 2.21 × 10⁻³) and validation (p = 3.63 × 10⁻⁶) cohorts, with consistent direction. The probe annotation has been clarified, and the Figure 3 / Table S5 discrepancy has been reconciled.
