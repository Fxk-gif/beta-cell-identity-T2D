# Response to Reviewer 2, Comment 7 — Xin identity-score verification + causal-language softening

## Reviewer comment (verbatim)

> (7) GSEA and protein-interaction analyses establish associations but cannot demonstrate an HNF1A-PDX1 regulatory cascade or transcriptionally driven β-cell failure. Notably, the full identity score was nonsignificant in the Xin single-cell dataset. The proposed "maturation-exhaustion" phenotype remains insufficiently distinguished from established β-cell dysfunction or dedifferentiation models.

## Findings

All three points raised by the reviewer are valid and are addressed as follows.

### 1. GSEA/PPI establish association, not a regulatory cascade (causal language softened)

The reviewer is correct. GSEA and PPI are associative analyses and cannot demonstrate an HNF1A-PDX1 regulatory cascade or transcriptionally driven β-cell failure. Causal wording in the original manuscript ("cascade", "transcriptionally driven", "regulatory lesion") has been revised to associative wording (consistent with Reviewer 1, Comment 7).

### 2. The full identity score in Xin is indeed nonsignificant

The full identity score (mean of per-gene z-scores across 16 identity genes) was recomputed in the Xin dataset (GSE81608, 18 donors, 503 β cells):

| Analysis level | T2D mean | ND mean | p value |
|---|---|---|---|
| Cell level (descriptive) | −0.0032 (n=296) | +0.0045 (n=207) | 0.9124 |
| **Donor level (primary inference)** | −0.0051 (n=6) | +0.0072 (n=12) | **0.8201** |

- Cohen's d (donor level) = −0.092, a negligible effect size
- **Conclusion:** the full identity score is indeed **nonsignificant** between T2D and ND; the reviewer's assessment is correct.

**Individual genes are nevertheless significant at the donor level** (with mixed directions that largely cancel in the aggregate score):

| Gene | Donor-level p value | Direction |
|---|---|---|
| GCK | 0.003 | up |
| PCSK2 | 0.003 | up |
| SLC2A2 | 0.018 | **down** |

> Note: the donor-level p values above are computed with the exact
> Mann–Whitney U test on per-donor mean raw RPKM expression (GCK/PCSK2 p=0.003,
> SLC2A2 p=0.018). Under the asymptotic normal approximation the values become
> GCK/PCSK2 p≈0.005 and SLC2A2 p≈0.019; the directions and significance
> conclusions are unchanged, with only minor numerical differences due to the
> approximation. The conclusions — SLC2A2 downregulated, GCK/PCSK2 upregulated —
> hold under both approaches.

This is consistent with the earlier finding that the single-cell signal is concentrated in individual genes (most consistently the SLC2A2 downregulation) rather than a global shift of the full identity score, and with the Comment 3 observation that concordance is driven by SLC2A2.

### 3. "Maturation–exhaustion" insufficiently distinguished from dedifferentiation models (conceptual boundary clarified)

The "maturation–exhaustion" phenotype overlaps substantially with established β-cell dysfunction / dedifferentiation models, and the data do not support it as a mechanistically independent entity. The manuscript now states that the findings should be interpreted as "consistent with established dedifferentiation/dysfunction frameworks" rather than as a separable new entity.

## Impact on the manuscript's core claims

**Moderate-to-substantial, in the direction of further narrowing the argument:**

1. The full identity score being nonsignificant at the single-cell level directly weakens the claim of a coordinated decline of the identity program at single-cell resolution; the single-cell evidence must be downgraded from "significant decline of the full score" to "downregulation of individual genes (SLC2A2)".
2. The causal-language softening (GSEA/PPI cannot demonstrate a cascade) is consistent with Reviewer 1, Comment 7, and is a necessary narrowing of the manuscript's wording.
3. The "maturation–exhaustion" concept is downgraded from an "independent mechanism" to "an expression pattern consistent with dedifferentiation", echoing the withdrawal of the "mature β-cell exhaustion" claim in Comment 4.
4. **What is retained:** the SLC2A2 single-gene downregulation (donor-level p=0.018; both bulk cohorts p<1e-5) and the significant bulk enrichment in the discovery cohort. This is the strongest and the only evidence that holds simultaneously across both bulk and single-cell analyses.

## Usage

```
python xin_identity_score_check.py
```

Dependencies: numpy, scipy. Data are located automatically from the project root under `sc_data/`.
