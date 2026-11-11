# Donor-level single-cell reanalysis — Reviewer 2, Comment 1

## Reviewer's comment

> "Individual cells were treated as independent biological replicates without
> accounting for donor-level variability. Critically, all 495 T2DM β-cells in
> GSE84133 originated from a single donor, rendering the extremely small
> reported P values potentially misleading. Donor-level pseudobulk or
> mixed-effects analyses and validation in additional diabetic donors are
> essential."

## Response summary

The reviewer is correct. We re-examined the donor structure of both single-cell
datasets and re-performed the differential-expression analysis at the donor
level, so that cells from the same donor are no longer treated as independent
observations. The key findings are:

1. **Baron et al. (GSE84133)**: the T2D group is a single donor (human4,
   495 β-cells). The cell-level P value reported in the original submission is
   therefore pseudoreplicated and has been **withdrawn**. For this dataset we
   now report only a directional log2 fold change (T2D vs. the mean of the three
   non-diabetic donors).

2. **Xin et al. (GSE81608)**: this dataset provides the additional diabetic
   donors requested (6 T2D vs. 12 non-diabetic donors). A donor-level two-tailed
   Wilcoxon rank-sum test with Benjamini–Hochberg correction was performed on
   per-donor mean expression.

## Files

| File | Description |
|---|---|
| `donor_level_reanalysis.py` | Python script performing both analyses |
| `baron_donor_structure.txt` | Beta-cell count per donor in the Baron dataset |
| `baron_donor_level_log2fc.csv` | Donor-level directional log2FC (Baron) |
| `xin_donor_level_results.csv` | Donor-level Wilcoxon + BH results (Xin) |

## How to run

```bash
# Set the data directory (contains the single-cell files) and run:
export T2D_SC_DATA=/path/to/sc_data
python donor_level_reanalysis.py
```

## Results

### Baron et al. (GSE84133) — donor structure

| Donor | Condition | Beta-cells |
|---|---|---|
| human1 | non-diabetic | 872 |
| human2 | non-diabetic | 371 |
| human3 | non-diabetic | 787 |
| human4 | type 2 diabetic | 495 |

The T2D group contains a single donor. No inferential test is valid; only a
directional log2 fold change is reported. SLC2A2 shows a marked directional
decrease (log2FC = −2.35).

### Xin et al. (GSE81608) — donor-level Wilcoxon + BH correction (18 donors)

| Gene | log2FC | p (donor-level) | FDR |
|---|---|---|---|
| SLC2A2 | −2.713 | 1.82×10⁻² | 0.097 |
| GCK | +1.334 | 3.23×10⁻³ | 0.026 * |
| PCSK2 | +0.651 | 3.23×10⁻³ | 0.026 * |

The remaining identity genes were not significant after correction. SLC2A2
retains a consistent directional decrease across both single-cell datasets and
across the two bulk cohorts (see the bulk confounder-adjustment analysis), but
is not independently significant at donor level in the single-cell data after
multiple-testing correction.
