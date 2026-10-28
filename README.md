# Analysis code for "Coordinated erosion of the β-cell identity program across independent cohorts defines a transcriptionally driven maturation-exhaustion phenotype in type 2 diabetic islets"

This repository contains the complete analysis code for the manuscript (Frontiers in Endocrinology, Manuscript ID 2015413). All statistical results reported in the paper and in the point-by-point response to reviewers can be reproduced by running these scripts.

## Authors

Xiankun Fu¹, Xiangyi Chen², Junling Jiang², Yong Xie²

¹Fujian University of Chinese Medicine, Fuzhou, Fujian, China  
²The Second Affiliated Hospital of Guangxi University of Chinese Medicine, Nanning, Guangxi, China

## Data sources

All datasets are publicly available from GEO:

| Accession | Description | Platform | Donors / Samples |
|---|---|---|---|
| GSE76895 | Human pancreatic islets, partial pancreatectomy (PPP) cohort | GPL570 (Affymetrix HG-U133 Plus 2.0) | 103 samples (36 T2D, 32 ND, 15 IGT, 20 T3cD) |
| GSE76894 | Human pancreatic islets, organ-donor (OD) cohort | GPL570 | 103 samples (19 T2D, 84 ND) |
| GSE81608 | Single-cell RNA-seq (Smart-seq2), Xin et al. | Smart-seq2 | 18 donors (12 ND, 6 T2D) |
| GSE84133 | Single-cell RNA-seq (inDrop), Baron et al. | inDrop | 4 donors (3 ND, 1 T2D) |

The two bulk datasets (GSE76895, GSE76894) originate from the IMIDIA biobank but derive from **distinct, non-overlapping donor groups** (partial-pancreatectomy vs. organ-donor), as stated in the manuscript.

## Environment

- **R**: version 4.6.1. Required packages: `ggplot2`, `pROC`, `tidyr`, `dplyr` (available from CRAN).
- **Python**: version 3.13. Required packages: `numpy`, `pandas`, `scipy` (available from PyPI).

## Script descriptions

### Main bulk analyses

| Script | Purpose |
|---|---|
| `analysis.R` | Defines the 19-gene β-cell identity score and the neuroendocrine/dedifferentiation scores; runs the main bulk differential analysis on GSE76895 (training) |
| `validation.R` | Independent cross-validation of the identity score in GSE76894 |
| `consistency.R` | Cross-dataset concordance of the identity score (core evidence figure) |

### Single-cell analyses (donor-level)

| Script | Purpose |
|---|---|
| `sc_validation.R` | Single-cell cross-validation in Xin et al. (GSE81608) |
| `sc_extract.py` | Extracts β-cell expression profiles from single-cell data |
| `sc_stats.py` | Single-cell β-cell statistics |
| `baron_extract.py` | Extracts β-cell profiles from Baron et al. (GSE84133) |

### Gene set enrichment (GSEA)

| Script | Purpose |
|---|---|
| `gsea.py` | Original GSEA implementation (hand-written enrichment score + 2,000 permutations) |
| `gsea_overlap_corrected.py` | **Overlap-corrected GSEA** — removes genes shared between curated TF target sets and the identity score, to address reviewer circularity concern (R1-3). Includes a sanity check that reproduces the original ES values exactly. |
| `gsea_plot.R` | GSEA enrichment plots |

### Confounder adjustment (reviewer R1-5)

| Script | Purpose |
|---|---|
| `confounder_adjustment.py` | Multiple linear regression `identity score ~ group + age + sex + BMI + batch` using sample-level covariates parsed from the GEO series matrices; reports adjusted and unadjusted p-values |

### Diagnostic models and networks

| Script | Purpose |
|---|---|
| `ml_identity_model.py` | β-cell identity diagnostic models (logistic regression / SVM / random forest) with ROC, calibration, and nomogram |
| `ml_identity_consistency.py` | Consistency analysis of the coordinated down-regulation of the identity network |
| `tf_network_figure.py` | Transcription-factor regulatory network (PDX1/HNF1A) |
| `ppi_network_figure.py` | STRING protein–protein interaction network |

### Supplementary analyses

| Script | Purpose |
|---|---|
| `reviewer_supplement_analysis.py` | Supplementary analyses: cell-autonomous core score (SLC2A2+IAPP+INS) and progenitor/disallowed score reversal |

### Data preprocessing

| Script | Purpose |
|---|---|
| `collapse_genes.py` | Collapses microarray probes to gene symbols |
| `build_r3_inputs.py` | Builds the analysis input files (`r3_inputs/`) |

## Reproduction notes

1. Download the four GEO datasets listed above.
2. Run `collapse_genes.py` and `build_r3_inputs.py` to prepare the input matrices.
3. Run `analysis.R` for the main bulk analysis, `validation.R` and `consistency.R` for validation and concordance.
4. Run `gsea.py` then `gsea_overlap_corrected.py` for GSEA (the latter reproduces the original ES values before removing shared genes).
5. Run `confounder_adjustment.py` for the confounder-adjusted results.

All p-values in the overlap-corrected GSEA and the confounder adjustment were computed with these scripts and match the values reported in the response to reviewers.

## License

This code is provided for reproducibility of the associated manuscript. Please cite the manuscript if you use it.
