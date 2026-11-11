# Response to Reviewer 2, Comment 6 — Gene-set reassessment

## Reviewer comment (verbatim)

> (6) The curated gene sets contain redundant identifiers (GLUT2/SLC2A2) and genes reportedly lacking microarray coverage. Moreover, Figure 9 incorrectly identifies SLC2A2 as the sole shared PDX1/HNF1A target, although INS and GCK are also shared. Gene-set definitions, effective gene counts, and enrichment statistics require careful reassessment.

## Findings

The three issues are addressed as follows.

### 1. Redundant identifiers (corrected)

| Gene set | Original list | After de-duplication | Removed redundant entry |
|---|---|---|---|
| PDX1 target genes | 10 | 9 | GLUT2 (= historical name for SLC2A2) |
| HNF1A target genes | 12 | 12 | none |
| β functional genes | 16 | 15 | PC1 (= historical name for PCSK1) |

### 2. Shared targets (Figure 9 corrected)

After de-duplication, PDX1 target genes ∩ HNF1A target genes = **{SLC2A2, INS, GCK}** (three genes), not SLC2A2 alone.

- PDX1 target genes (9 after de-duplication): G6PC2, GCK, IAPP, INS, MAFA, NKX6-1, PAX4, PC, SLC2A2
- HNF1A target genes (12): ABCC8, FOXA3, GCK, HNF1B, HNF4A, HNF4G, INS, KCNJ11, PDX1, PKLR, SLC2A2, SLC30A8
- Intersection: **GCK, INS, SLC2A2**

Figure 9 has been corrected to mark the three shared targets.

### 3. Genes lacking microarray coverage (probe-coverage check)

Probe-coverage check on the GPL570 (HG-U133 Plus 2.0) platform:

- **MAFA: 0 probes** (no coverage; excluded from the effective gene count)
- GLUT2: 0 probes (an alias of SLC2A2, merged into SLC2A2 after de-duplication; SLC2A2 has 2 probes)
- All other genes have ≥1 probe

### 4. Effective gene counts (after de-duplication + coverage filtering)

| Gene set | Listed | Effective | Excluded |
|---|---|---|---|
| PDX1_targets | 9 | 8 | MAFA |
| HNF1A_targets | 12 | 12 | none |
| beta_functional | 15 | 15 | none |

### 5. GSEA enrichment statistics (correct log2FC, recomputed after de-duplication + coverage filtering)

**Methodological note:** the expression matrix is in RMA-normalized log2 space; fold change is computed as the difference in group means (m_T2D − m_ND), not log2(ratio).

| Gene set | Discovery cohort GSE76895 (ES, p) | Validation cohort GSE76894 (ES, p) |
|---|---|---|
| PDX1_targets | −0.948, p<0.001 | −0.892, p=0.074 |
| HNF1A_targets | −0.808, p=0.015 | −0.687, p=0.907 |
| beta_functional | −0.815, p=0.005 | −0.902, p=0.011 |

## Impact on the manuscript's core claims

**Moderate, with a directionally contained effect:**

1. The de-duplication and shared-target correction are purely corrective: they do not change the direction of the conclusions and make Figure 9 and the gene-set definitions more accurate.
2. In the **discovery cohort**, all three gene sets retain significant downregulation enrichment (PDX1 p<0.001, HNF1A p=0.015, functional p=0.005), so the core claim that PDX1/HNF1A target genes are enriched among downregulated genes remains valid in the discovery cohort.
3. In the **validation cohort**, PDX1 (p=0.074) and HNF1A (p=0.907) enrichment are not significant — consistent with the Comment 3 finding that cross-cohort concordance is driven by SLC2A2. Only the beta_functional gene set (including SLC2A2, GCK, etc.) remains significant in the validation cohort.
4. Net effect: the "PDX1/HNF1A cascade" claim is further narrowed to the discovery cohort plus the SLC2A2/functional-gene level, consistent with the downgraded wording used elsewhere in the response.

## Usage

```
python gene_set_reassessment.py
```

Dependencies: numpy, scipy. Data are located automatically from the project root under `r3_inputs/` and `geo_data/processed/`.
