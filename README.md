# beta-cell-identity-T2D — Analysis code (revision)

Analysis code for the manuscript *"Reduced β-cell identity across separate datasets defines an expression pattern of dedifferentiation in type 2 diabetic islets"*.

This repository contains the scripts used for the revised analyses, organised by the reviewer comments they address. All analyses run on publicly available data (GEO), in a Python 3 environment. Input data are located via environment variables.

## Repository structure

```
.
├── reviewer1/   # Reviewer 1, Comments 1–8 (code modules)
├── reviewer2/   # Reviewer 2, Comments 1–8 (code modules)
└── figures/     # Scripts that regenerate the revised figures
```

### `reviewer1/`

| Directory | Reviewer 1 comment | Analysis |
|---|---|---|
| `01_donor_level_single_cell/` | #1 pseudoreplication | Donor-level single-cell re-analysis (Baron GSE84133, Xin GSE81608) |
| `02_donor_overlap_check/` | #2 "independent cohorts" | Three-layer donor-overlap check (GSE76895 vs GSE76894) |
| `03_gsea_circularity_correction/` | #3 GSEA circularity | Overlap-corrected / de-circularized GSEA |
| `04_cell_fraction_deconvolution/` | #4 cell fraction | Per-donor cell-type fraction from single-cell annotations |
| `05_confounder_adjustment/` | #5 confounders | Identity score ~ group + age + sex + BMI + batch (OLS) |
| `06_diagnostic_overfitting/` | #6 overfitting | CV, DeLong test, bootstrap CI, class-imbalance metrics |
| `07_progenitor_trace_probes/` | #8 trace-level probes | Bulk progenitor/disallowed score at full precision |

### `reviewer2/`

| Directory | Reviewer 2 comment | Analysis |
|---|---|---|
| `01_donor_level_reanalysis/` | #1 donor-level statistics | Baron + Xin donor-level pseudobulk re-analysis |
| `02_sensitivity_analysis/` | #2 sensitivity | Sensitivity analysis of the identity-score results |
| `03_cross_cohort_correlation/` | #3 cross-cohort | Leave-one-out and rank-based cross-cohort correlation |
| `04_maturation_state_quantification/` | #4 maturation state | Per-cell maturation score and mature-β proportion |
| `05_probe_annotation_uniform_panel/` | #5 probe annotation | Uniform 15-gene panel across cohorts |
| `06_gene_set_reassessment/` | #6 gene-set redundancy | Gene-set reassessment (redundant IDs, coverage) |
| `07_identity_score_check/` | #7 Xin identity score | Full identity-score re-computation in Xin |
| `08_roc_uncertainty_stability/` | #8 ROC uncertainty | Bootstrap CI, DeLong, repeated CV, imbalance metrics |

## Data sources

| Data | GEO accession |
|---|---|
| Bulk microarray, discovery | GSE76895 (IMIDIA, GPL570) |
| Bulk microarray, validation | GSE76894 (IMIDIA, GPL570) |
| Single-cell, Xin et al. 2016 | GSE81608 |
| Single-cell, Baron et al. 2016 | GSE84133 |

The scripts expect the following input directories (paths resolved via environment variables or by walking upward from the script location):

- `r3_inputs/` — bulk expression matrices (`GSE7689*_expr.csv`, `_groups.csv`), TF target sets, volcano input
- `geo_data/` — GEO series-matrix files (`GSE7689*_series_matrix.txt.gz`)
- `sc_data/` — single-cell data (`xin_beta_expr.csv`, `xin_beta_annot.csv`, Baron count files)

## Running the analyses

```bash
# Point the environment variables at the workspace root (optional; the scripts
# also locate the directories by walking upward from their own location)
export T2D_DATA=/path/to/workspace           # directory containing r3_inputs/, geo_data/, sc_data/
export T2D_SC_DATA=/path/to/workspace/sc_data
export T2D_R3_INPUTS=/path/to/workspace/r3_inputs

# Reviewer 1
python reviewer1/01_donor_level_single_cell/xin_donor_level.py
python reviewer1/01_donor_level_single_cell/baron_donor_level.py
python reviewer1/02_donor_overlap_check/donor_overlap_check.py
python reviewer1/03_gsea_circularity_correction/overlap_corrected_gsea.py
python reviewer1/04_cell_fraction_deconvolution/cell_fraction_analysis.py
python reviewer1/05_confounder_adjustment/confounder_adjustment.py
python reviewer1/06_diagnostic_overfitting/overfitting_assessment.py
python reviewer1/07_progenitor_trace_probes/progenitor_score_bulk.py

# Reviewer 2
python reviewer2/01_donor_level_reanalysis/donor_level_reanalysis.py
python reviewer2/02_sensitivity_analysis/sensitivity_analysis.py
python reviewer2/03_cross_cohort_correlation/leave_one_out_correlation.py
python reviewer2/04_maturation_state_quantification/maturation_state_quantification.py
python reviewer2/05_probe_annotation_uniform_panel/uniform_panel_reanalysis.py
python reviewer2/06_gene_set_reassessment/gene_set_reassessment.py
python reviewer2/07_identity_score_check/xin_identity_score_check.py
python reviewer2/08_roc_uncertainty_stability/roc_uncertainty_and_model_stability.py
```

## Figure scripts

The `figures/` directory contains the scripts that regenerate the six figures
revised during peer review (Figure 2, 3, 6, 7 and Supplementary Figures 1, 3).
See `figures/README.md` for the figure-to-script mapping and the R/Python
requirements. These scripts locate the analysis workspace via the
`T2D_WORKSPACE` environment variable (or by walking upward to a directory
containing `r3_inputs/`).

## Dependencies

- Python 3
- numpy, pandas, scipy, scikit-learn
- R 4.x (for the figure scripts): ggplot2, pROC, gridExtra

Each module directory contains a `README.md` with the corresponding reviewer comment and a description of the analysis, together with the output files produced by running the script.
