# Reviewer 1, Comment #8 — Trace-level probes and the progenitor/disallowed score

## Reviewer comment

> The argument against dedifferentiation rests on trace-level probes. The
> progenitor/disallowed score uses eight genes the authors themselves describe
> as expressed at trace levels in bulk islets. A p value of 1.4×10⁻¹¹ on
> near-floor microarray signal is weak evidence against an active
> dedifferentiation program.

## Files in this folder

| File | Description |
|---|---|
| `progenitor_score_bulk.py` | Recomputation of the bulk progenitor/disallowed score (8-gene panel, exact Wilcoxon) |
| `progenitor_score_output.txt` | Complete output from running the script |
| `progenitor_score_results.csv` | Score-level and single-gene results table |

## What the script does

The reviewer's concern is that a highly significant p-value computed on genes
expressed at trace levels may reflect the variance structure of the microarray
platform rather than biological signal. To make the basis of the claim
inspectable, this script recomputes the bulk progenitor/disallowed score at
full precision and reports each gene's group-mean expression and exact p-value.

- **Score**: mean of per-gene z-scores across samples, for the 8-gene panel
  `FOXE1, MYCN, CHL1, OPRD1, GRIN2D, ALDH1A3, HES1, SOX9`.
- **Test**: Wilcoxon rank-sum (Mann-Whitney U), exact method (matching
  `wilcox.test(..., exact = TRUE)` in R).
- **Platform note**: `FOXE1` has no probe on the GPL570 platform of the
  validation cohort (GSE76894), so the validation score is computed from the
  remaining 7 genes; this is printed in the output.

## Usage

```bash
# Point T2D_DATA at the project directory that contains r3_inputs/
export T2D_DATA=/path/to/project
python progenitor_score_bulk.py
```

Dependencies: Python 3 + numpy + scipy.

## Results

| Cohort | Genes used | Score (ND vs T2D) | Exact p |
|---|---|---|---|
| GSE76895 (discovery) | 8 | +0.321 vs −0.286 | **1.37×10⁻¹¹** |
| GSE76894 (validation) | 7 (no FOXE1) | +0.052 vs −0.229 | **2.45×10⁻⁴** |

The score is lower in T2D than in non-diabetic islets in both cohorts, a
direction opposite to what an active dedifferentiation program would predict.
However, several of the eight genes are expressed at low (trace) levels in bulk
islets (e.g. group means near the floor of the RMA scale), so this reduction is
read as reflecting the same depletion of mature β-cells rather than as strong
evidence against dedifferentiation.

## Response to the manuscript

In line with the reviewer's point, the manuscript no longer presents the
progenitor/disallowed score as strong evidence against an active
dedifferentiation program. The Discussion instead states that the bulk-level
data "provide no clear evidence of an active dedifferentiation program
detectable at bulk resolution", and the trace-level caveat is reported
alongside the broader, converging evidence (donor-level single-cell re-analysis
and the direction of the score change). The reported p-values are standardised
to full precision, matching Supplementary Table S8B.
