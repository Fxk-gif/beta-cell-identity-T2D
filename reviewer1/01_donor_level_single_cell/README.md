# Donor-level single-cell reanalysis

This directory contains the analysis scripts used to correct the
pseudoreplication issue in the single-cell analyses.

## Background

Reviewer 1, Comment 1 noted that the original manuscript ran Mann-Whitney tests
treating individual cells as independent samples, which is pseudoreplication
because multiple beta cells from the same donor are not independent. In
particular, the T2DM group of the Baron dataset (GSE84133) consists of a single
donor, yet FDR = 2.4e-34 was reported and Figure 5 was titled "Independent
replication". These scripts aggregate the unit of analysis from individual cells
to donors (donor-level pseudobulk).

## Files

| File | Description |
|---|---|
| `xin_donor_level.py` | Xin (GSE81608) donor-level reanalysis: 18 donors (12 ND + 6 T2D), Wilcoxon test with BH correction |
| `baron_donor_level.py` | Baron (GSE84133) donor-level reanalysis: T2D has only 1 donor, so only a directional log2FC is reported |

## Data

Data must be downloaded from GEO and placed in the `sc_data/` directory. The
scripts locate this directory via the `T2D_SC_DATA` environment variable, or by
walking upward from the script location to find the project root that contains
`sc_data/`.

- **Xin (GSE81608)**: `GSE81608_human_islets_rpkm.txt.gz` (RPKM expression
  matrix). Beta cells must first be extracted to produce `xin_beta_expr.csv`
  (genes x cells) and `xin_beta_annot.csv` (cell -> donor and disease status).
- **Baron (GSE84133)**: `GSM2230757`-`GSM2230760` (human1-4) UMIFM counts
  (human1-3 non-diabetic, human4 T2D).

## Usage

```bash
# Optional: set the data directory (otherwise it is located by walking upward)
export T2D_SC_DATA=/path/to/sc_data

python xin_donor_level.py      # Xin: donor-level Wilcoxon + BH FDR
python baron_donor_level.py    # Baron: donor-level directional log2FC
```

Dependencies: Python 3 + numpy (no other third-party packages).

## Results

### Xin (GSE81608, 18 donors: 12 ND + 6 T2D)

Donor-level Wilcoxon rank-sum test with Benjamini-Hochberg correction:

| Gene | log2FC | p (donor-level) | FDR | Result |
|---|---|---|---|---|
| SLC2A2 | -2.71 | 0.018 | 0.103 | not significant (cell-level pseudoreplication corrected) |
| IAPP | -0.44 | 0.68 | 0.83 | not significant |
| GCK | +1.33 | 0.003 | 0.027 | significantly up-regulated |
| PCSK2 | +0.65 | 0.003 | 0.027 | significantly up-regulated |

The apparent down-regulation of SLC2A2 and IAPP in the cell-level test does not
replicate at the donor level, consistent with the pseudoreplication concern.
At the donor level, GCK and PCSK2 are instead significantly up-regulated.

### Baron (GSE84133, T2D has only 1 donor)

| Gene | ND mean counts | T2D counts | log2FC (directional only) |
|---|---|---|---|
| SLC2A2 | 32.7 | 6.0 | -2.35 |
| IAPP | 134853 | 50749 | -1.41 |
| INS | 730437 | 384175 | -0.93 |

Because the T2D group contains only 1 donor (human4, 495 beta cells), no
p-value or FDR is reported; only a directional log2FC is given, consistent with
the principle that a fold change is the only defensible statistic when n = 1.

## Notes

- The public Xin data are RPKM (not raw counts), so an edgeR negative binomial
  model does not apply; a donor-level non-parametric Wilcoxon test is used
  (unit = donor, n = 18).
- The Wilcoxon test uses the exact distribution
  (`scipy.stats.mannwhitneyu(method='exact')`), with a normal approximation and
  tie correction as a fallback when scipy is unavailable. The BH FDR correction
  is implemented explicitly.
