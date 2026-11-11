# Figure scripts for the revised manuscript

This folder contains the scripts that regenerate the six figures that were
revised during peer review (the remaining figures are unchanged from the
initial submission).

| Script | Figure | Language |
|---|---|---|
| `figure2_identity_score.R` | Figure 2 — reduced beta-cell identity score (box plots) | R |
| `figure3_roc.R` | Figure 3 — diagnostic ROC curves of identity genes | R |
| `figure6_cross_dataset.py` | Figure 6 — cross-dataset concordance of single-cell changes | Python |
| `figure7_gsea.py` | Figure 7 — GSEA of TF target gene sets (original + overlap-corrected) | Python |
| `suppfig1_heatmap.py` | Supplementary Figure 1 — group-mean expression heatmap | Python |
| `suppfig3_ppi.py` | Supplementary Figure 3 — STRING PPI network | Python |

## Data location

The scripts read their inputs from the analysis workspace, not from this
folder. The workspace is resolved in this order:

1. the `T2D_WORKSPACE` environment variable, if set; otherwise
2. by walking upward from the script directory to find a folder that contains
   `r3_inputs/`.

The workspace contains the following input files:

- `r3_inputs/GSE76895_expr.csv`, `r3_inputs/GSE76895_groups.csv`
- `r3_inputs/GSE76894_expr.csv`, `r3_inputs/GSE76894_groups.csv`
- `r3_inputs/identity_gene_sets.csv`
- `r3_inputs/volcano_input.csv`, `r3_inputs/tf_target_sets.csv`
- `sc_data/xin_beta_expr.csv`, `sc_data/xin_beta_annot.csv`
- `sc_data/baron_donor_level_counts.csv`
- `ml_output/ppi_string_result.xlsx`
- `ml_output/identity_log2fc_both.csv`

The raw single-cell counts and bulk matrices were obtained from the Gene
Expression Omnibus (GSE81608, GSE84133, GSE76895, GSE76894); the STRING
interaction table was downloaded from the STRING database (v12).

## Requirements

- R 4.6.1 with packages `ggplot2`, `pROC`, `gridExtra`
- Python 3.13 with `numpy`, `scipy`, `pandas`, `matplotlib`

## Running

```bash
# set the workspace (optional if r3_inputs/ is reachable by walking upward)
export T2D_WORKSPACE=/path/to/workspace

Rscript figure2_identity_score.R
Rscript figure3_roc.R

python figure6_cross_dataset.py
python figure7_gsea.py
python suppfig1_heatmap.py
python suppfig3_ppi.py
```

Output figures are written to `figures/` inside the workspace as
`Figure2.TIF`, `Figure3.TIF`, `Figure6.TIF`, `Figure7.TIF`,
`SupplementaryFigure1.TIF`, and `SupplementaryFigure3.TIF`.

## Statistical notes

- Figure 2 group comparison uses the exact Wilcoxon rank-sum test
  (`wilcox.test(..., exact = TRUE)`), matching the values reported in the
  manuscript (discovery p = 2.21e-3, validation p = 3.63e-6).
- Figure 6 Xin statistics are donor-level Mann-Whitney U tests with
  Benjamini-Hochberg FDR correction; the Baron cohort has a single T2D donor,
  so its changes are reported descriptively only.
- Figure 7 enrichment scores are computed with the weighted running-sum
  statistic (seed = 42, 2000 permutations), identical to
  `overlap_corrected_gsea.py`.
