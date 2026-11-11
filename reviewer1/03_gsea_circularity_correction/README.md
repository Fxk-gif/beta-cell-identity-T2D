# Reviewer 1, Comment #3 — GSEA circularity verification and de-circularized re-run

> Corresponding reviewer comment (Reviewer 1, Comment #3):
> *Partial circularity between the identity score and the GSEA. The curated
> PDX1 and HNF1A target sets share members with the 16-gene score (SLC2A2, GCK,
> IAPP, INS, NKX6-1, PDX1). The cascade is also in part a curation choice: PDX1
> was placed in the HNF1A target list, after which its downstream position is
> presented as a structural finding (Figure 9). Separately, the hypothesis-free
> analysis returned KLF10 rather than PDX1 or HNF1A, which reads more as
> non-replication of the curated result than as the confirmation it is framed
> to be.*

---

## Files in this folder

| File | Description |
|---|---|
| `overlap_corrected_gsea.py` | Circularity quantification + de-circularized GSEA re-run script (based on the actual data) |
| `overlap_corrected_gsea_results.csv` | Summary table of the de-circularized results |
| `overlap_corrected_gsea_output.txt` | Complete output from running the de-circularized script |

---

## What the script does (all based on the actual data)

1. **Quantify circularity** — compute the genes shared between the curated
   target gene sets and the 16-gene identity score
2. **De-circularized GSEA** — after removing the identity-score genes from each
   target set, re-run using the exact same algorithm as the original `gsea.py`
   (seed=42, 2000 permutations) and report the change in significance
3. **Verify curation facts** — check whether PDX1 appears in the HNF1A target
   list (and vice versa)

---

## Results

### (1) Circularity quantification — the reviewer's concern is confirmed

| Target gene set | Size | Shared with 16-gene score | Shared genes |
|---|---|---|---|
| PDX1_targets | 9 | 5 | GCK, IAPP, INS, NKX6-1, SLC2A2 |
| HNF1A_targets | 12 | 7 | ABCC8, GCK, HNF4A, INS, PDX1, SLC2A2, SLC30A8 |
| beta_identity_network | 18 | 10 | GCK, HNF1A, HNF4A, IAPP, INS, NEUROD1, NKX6-1, PAX6, PDX1, SLC2A2 |
| beta_functional_genes | 15 | 9 | ABCC8, GCK, GLP1R, IAPP, INS, PCSK1, PCSK2, SLC2A2, SLC30A8 |

> All six shared genes named by the reviewer (SLC2A2, GCK, IAPP, INS, NKX6-1,
> PDX1) are confirmed to be present.

### (2) De-circularized GSEA — the enrichment signal largely disappears (important)

| Target set | Original ES (p) | Genes retained after de-circularization | De-circularized ES (p) | Conclusion |
|---|---|---|---|---|
| PDX1_targets | −0.948 (p<0.001) | 2 genes | — | **cannot be tested (<3 genes)** |
| HNF1A_targets | −0.808 (p=0.015) | 4 | +0.378 (p=0.966) | **not significant, direction reversed** |
| beta_identity_network | −0.780 (p=0.012) | 7 | −0.276 (p=0.997) | **not significant** |
| beta_functional_genes | −0.815 (p=0.005) | 4 | −0.901 (p=0.034) | still significant (downregulated) |

> **Core conclusion: after removing the identity-score genes, the enrichment
> signal of the PDX1 and HNF1A target sets disappears (the PDX1 set cannot be
> tested; the HNF1A set is no longer significant; only the β-cell functional
> gene module remains significantly downregulated). The reviewer's circularity
> concern is valid.**

### (3) Cascade curation facts — the reviewer is correct

- `PDX1 ∈ HNF1A_targets`: **Yes** (PDX1 is indeed placed in the HNF1A target list)
- `HNF1A ∈ PDX1_targets`: No (unidirectional, not mutual)

---

## Regarding KLF10 (hypothesis-free analysis) — positioned as a curation-independent finding

The reviewer noted that "hypothesis-free analysis returned KLF10 rather than
PDX1 or HNF1A". The original manuscript reported: fgsea + hTFtarget database,
KLF10 NES = −1.977, BH-adjusted p = 5.7×10⁻⁷, ranked first (PDX1/HNF1A did not
rank near the top).

**Positioning**: the KLF10 result comes from the hTFtarget database; the
script for this hypothesis-free analysis is provided here. KLF10 is treated as a
curation-independent hypothesis-free finding, used only to indicate a
stress-signalling (TGF-β/KLF10) component in the phenotype, and **does not serve
as independent validation of the PDX1/HNF1A cascade**.

Accordingly, the analysis script for this comment contains only the circularity
quantification and the de-circularized GSEA (`overlap_corrected_gsea.py`); it
does not re-implement the hTFtarget hypothesis-free analysis. This scope is
consistent with the response document and the final manuscript.

---

## Revisions to the manuscript argument (made in the response letter)

Based on the results above, the response letter makes the following changes:

1. Acknowledge the circularity and remove the statement that "GSEA serves as
   independent validation of the identity score";
2. Report explicitly that the enrichment signal of the PDX1/HNF1A target sets
   disappears after de-circularization;
3. Reframe Figure 9 as a "literature-based curation schematic" and remove the
   "structural finding" wording;
4. Re-describe the KLF10 hypothesis-free analysis as a "curation-independent
   finding" rather than a "confirmation", stating its hypothesis-free nature and
   not using it as independent support for the curated cascade.
