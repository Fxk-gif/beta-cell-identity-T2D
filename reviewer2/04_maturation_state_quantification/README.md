# Response to Reviewer 2, Comment 4 — Direct quantification of β-cell maturation state

## What this folder contains

| File | Description |
|---|---|
| `maturation_state_quantification.py` | Quantification of per-cell β-cell maturation score and donor-level (donor-adjusted) comparison, in Xin GSE81608. |
| `maturation_score_donor_level.csv` | Per-donor mean maturation score. |
| `mature_beta_proportion_donor_level.csv` | Per-donor mature β-cell proportion. |
| `README.md` | This document. |

## The reviewer's comment

> The absence of significant intracellular changes in PDX1 and NKX6-1 cannot establish depletion of mature β-cell subpopulations. This interpretation requires direct quantification of β-cell maturation states and their donor-adjusted proportions rather than inference from discordant bulk and single-cell expression patterns.

## What we did

We agree with the reviewer: the previous inference of mature β-cell depletion from the absence of PDX1/NKX6-1 changes was not justified. We therefore replaced it with a direct, donor-adjusted quantification.

**Maturation score** (per β-cell):

```
maturation score = mean(z-score of mature markers) − mean(z-score of immature/progenitor markers)
```

- Mature markers: MAFA, UCN3, INS, SLC2A2, GCK, SLC30A8, PCSK1 (SIX3 was not present in the Xin matrix)
- Immature/progenitor markers: MAFB, NEUROG3, CHL1, MYCN, OPRD1

z-scores are computed per gene across all β-cells. Genes with zero variance (e.g. NEUROG3, which is not expressed in mature β-cells) contribute 0.

**Donor-level (primary) analysis**: each donor's β-cells were aggregated to a mean maturation score, and T2D donors (n = 6) were compared with ND donors (n = 12) using a Mann-Whitney U test, treating the donor (not the cell) as the unit of analysis.

**Mature β-cell proportion**: each β-cell was classified as mature if its score exceeded the median ND score; the per-donor mature fraction was then compared between T2D and ND.

## Results (Xin, GSE81608)

### Donor-level maturation score

| Group | n | Mean score | SD |
|---|---|---|---|
| T2D donors | 6 | +0.0017 | 0.1226 |
| ND donors | 12 | −0.1037 | 0.2150 |

Mann-Whitney U (donor-level): **p = 0.3845** (not significant).

### Donor-level mature β-cell proportion

| Group | n | Mean mature fraction |
|---|---|---|
| T2D donors | 6 | 59.11% |
| ND donors | 12 | 48.16% |

Mann-Whitney U (donor-level): **p = 0.2908** (not significant).

### Individual maturation markers (cell-level, T2D vs ND)

| Marker | T2D mean | ND mean | log2FC | p |
|---|---|---|---|---|
| MAFA (mature) | 11.89 | 8.06 | +0.56 | 0.011 |
| MAFB (immature) | 47.09 | 53.32 | −0.18 | 0.719 |
| UCN3 (mature) | 3.70 | 2.51 | +0.56 | 0.315 |
| CHL1 (immature) | 6.86 | 11.39 | −0.73 | 0.0002 |

## Conclusion

At the donor level (the appropriate unit of analysis), there is **no significant reduction** in β-cell maturation state (p = 0.3845) or in the mature β-cell proportion (p = 0.2908) in T2D. If anything, the direction is toward *higher* maturity in T2D (MAFA upregulated, CHL1 downregulated). The term "depletion of mature β-cell subpopulations" is therefore not supported by these data and has been withdrawn from the manuscript.

Note: Baron GSE84133 contains only one T2D donor (vs three ND donors), so a donor-level statistical comparison is not possible for that dataset; it is reported descriptively only.
