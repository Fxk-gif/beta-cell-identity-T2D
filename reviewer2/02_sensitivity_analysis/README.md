# Response to Reviewer 2, Comment 2 — Sensitivity analysis

## Overview

| File | Description |
|---|---|
| `sensitivity_analysis.py` | Main sensitivity-analysis script (based on public GEO data) |
| `sensitivity_analysis_results.txt` | Output produced by running the script |

## Reviewer comment addressed

Reviewer 2, Comment 2:

> The bulk transcriptomic analyses do not adequately account for differences in β-cell abundance, islet purity, or clinical confounders. Furthermore, both bulk datasets originate from the same IMIDIA study and use the same microarray platform, limiting claims of cross-platform reproducibility. Cell-composition adjustment and appropriate sensitivity analyses are necessary.

## What this script does (four analyses)

1. **Glycemia adjustment.** Fasting glucose (GSE76895) or ICU glycemia (GSE76894) is added as a covariate to the OLS adjustment model to test whether the T2D effect on the identity score remains robust after further adjustment.

2. **β-cell abundance proxy adjustment.** The GEO metadata do not contain islet-purity or β-cell-proportion fields (a limitation noted by the reviewer), so a β-cell relative-abundance proxy is constructed from the mean expression of a set of islet endocrine-lineage marker genes (CHGA, CHGB, SCG2, SCG3, ISL1, PCSK1N) that are not part of the identity-score definition, and added as a covariate. The proxy is derived from the same bulk expression data and is an approximation rather than a formal deconvolution; it is labelled as a proxy throughout.

3. **Transcription-factor vs. effector gene sensitivity.** The 15-gene identity score (the 16-gene panel minus PDX1) is split into a 6-gene transcription-factor (TF) sub-score and a 9-gene effector/functional sub-score to test which class drives the T2D effect. If reduced β-cell proportion were the main driver, both classes should decline in parallel; a selectively stronger TF decline would instead support cell-autonomous erosion of the regulatory program.

4. **Same-origin / same-platform statement.** Confirms that the two cohorts belong to the same IMIDIA study and use the same GPL570 platform (as the reviewer noted), and records the corresponding downgraded wording in the output.

## Data sources

- GSE76895 (partial pancreatectomy, 68 samples: T2D = 36, ND = 32)
- GSE76894 (organ donor, 103 samples: T2D = 19, ND = 84)
- Both are from the IMIDIA study and use the Affymetrix GPL570 microarray platform.

## Usage

```bash
python sensitivity_analysis.py
```

The script locates the workspace directory containing `geo_data/` and `r3_inputs/` either through the `T2D_DATA` environment variable or by walking up from the script. Requires numpy.

## Results

| Cohort | Base-adjusted p | + glycemia p | + glycemia + β proxy p |
|---|---|---|---|
| GSE76895 (discovery) | 2.85e-1 | 3.55e-1 | 6.97e-3 |
| GSE76894 (validation) | 6.30e-3 | 3.46e-2 | 1.03e-1 |

The β-cell abundance proxy is significantly lower in T2D vs ND in GSE76894 (Mann-Whitney p = 3.67e-4) but not in GSE76895 (p = 1.73e-1).

**Main finding:** After adjustment for the β-cell abundance proxy, the identity-score difference in the validation cohort GSE76894 changes from significant (p = 6.3e-3) to nonsignificant (p = 1.03e-1), indicating that part of the bulk-level identity-score decline in the validation cohort is driven by β-cell abundance (compositional) confounding rather than a purely cell-autonomous effect. In the discovery cohort GSE76895, the identity score is nonsignificant after base adjustment (p = 2.85e-1) and after further glycemia adjustment (p = 3.55e-1), but the addition of the β-cell abundance proxy renders the association nominally significant (p = 6.97e-3), a pattern opposite to that in the validation cohort. The revised manuscript softens the corresponding conclusion accordingly and relies on the single-cell data as the primary source of cell-autonomous evidence.
