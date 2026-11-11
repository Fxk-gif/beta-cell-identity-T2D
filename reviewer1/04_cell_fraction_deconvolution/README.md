# Reviewer 1, Comment #4 — Cell-fraction analysis code

## Reviewer comment

> No deconvolution was performed, despite the term appearing in the Abstract. Compositional depletion of mature β-cells is inferred from the absence of a single-cell PDX1/NKX6-1 signal - absence of evidence, with dropout and low transcription-factor abundance as competing explanations. Segerstolpe et al. is cited as a composition reference, but no cell fractions are reported anywhere.

## Files in this folder

| File | Description |
|---|---|
| `cell_fraction_analysis.py` | Cell-type fraction calculation script (based on public data) |
| `cell_fraction_results.csv` | Counts and fractions of each cell type per donor |
| `cell_fraction_analysis_output.txt` | Complete output from running the script |

## Usage

```bash
# Set the data directory environment variable (pointing to the directory containing the Baron count files)
export T2D_SC_DATA=/path/to/sc_data
python cell_fraction_analysis.py
```

Dependencies: Python standard library only (`gzip`/`csv`/`numpy`; numpy is used
only for the mean/standard deviation and can be replaced).

## Data sources

- **Baron et al. (GSE84133)**: human islet inDrop single-cell dataset, 4 donors
  - `human1` / `human2` / `human3` = non-diabetic (ND)
  - `human4` = type 2 diabetes (T2D)
  - Column 3 (`assigned_cluster`) of each donor's `*_umifm_counts.csv.gz` holds the cell-type annotation
- **Xin et al. (GSE81608)**: purified β-cell dataset (all `subtype` = "beta"),
  **contains no other cell types and therefore cannot be used for fraction
  analysis** — stated as is.

## Results

**The cell-type annotations of the single-cell datasets do NOT support a
decrease in the β-cell (compositional) fraction.**

| Donor | Disease | Total cells | β of total cells | β of endocrine cells |
|---|---|---|---|---|
| human1 | ND | 1937 | 45.02% | 62.06% |
| human2 | ND | 1724 | 21.52% | 29.44% |
| human3 | ND | 3605 | 21.83% | 37.19% |
| **human4** | **T2D** | **1303** | **37.99%** | **52.44%** |

Key facts:
1. The β fraction of the T2D donor (human4) (37.99%) is **higher** than the
   mean of the ND donors (29.46%), opposite to the direction implied by
   "compositional depletion".
2. The β fraction varies substantially among the three ND donors (21.52%–45.02%;
   standard deviation 11 percentage points), far exceeding the ND–T2D
   difference, indicating that this variation is driven mainly by
   inter-donor/technical batch effects (islet isolation, cell sorting) rather
   than disease state.
3. The T2D group contains only 1 donor, so no statistical test is possible.

## Response to the manuscript

The core judgment in Comment #4 is **correct**: the manuscript inferred
"compositional depletion of mature β-cells" from the absence of a single-cell
PDX1/NKX6-1 signal — an absence-of-evidence argument — and the computed cell
fractions provide no evidence of a decreased β fraction.

The appropriate actions (implemented in the response document):
1. **Withdraw** the term "deconvolution" from the Abstract and main text (no
   computational deconvolution was performed), replacing it with an accurate
   description of the analysis actually performed;
2. **Withdraw/qualify** the compositional causal claim of "compositional
   depletion of mature β-cells";
3. Report the cell-fraction table above (as a new supplementary table) and
   state clearly that the single-cell data do not support a decreased β fraction;
4. Re-attribute the bulk decrease in PDX1/NKX6-1 to "dropout / low abundance /
   donor batch confounding cannot be excluded" rather than "reduced mature
   β-cells";
5. List dropout, low transcription-factor abundance, and donor batch as
   competing explanations in the Limitations section.
