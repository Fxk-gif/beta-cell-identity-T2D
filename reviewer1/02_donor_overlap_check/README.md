# Reviewer 1, Comment #2 — Donor overlap check code

> Corresponding reviewer comment (Reviewer 1, Comment #2):
> *"Independent cohorts" overstates the design. GSE76895 and GSE76894 are both
> IMIDIA biobank, both GPL570. Concordance between them cannot exclude biobank
> or platform artefact — which is precisely what the manuscript claims it
> establishes. Donor overlap should be checked and stated explicitly.*

---

## Files in this folder

| File | Description |
|---|---|
| `donor_overlap_check.py` | Donor overlap check script (Python standard library only, no third-party dependencies) |
| `donor_overlap_check_output.txt` | Complete output from running the script |

---

## What the script does

In response to the reviewer's request to "check and state explicitly" the donor
overlap between the two cohorts, the script performs a **three-layer overlap
check** on the public GEO metadata for GSE76895 (partial pancreatectomy [PPP]
cohort) and GSE76894 (organ donor [OD] cohort):

1. **Donor identifier cross-comparison** — extracts the donor code from each
   sample title:
   - GSE76895 donor code prefix `DP###` (pancreatectomy patients), numbering range 2–197
   - GSE76894 donor code prefix `PO###` (organ donors), numbering range 755–946
   - The two donor-code prefixes differ, and the numeric ranges are completely disjoint

2. **Numeric-ID cross-comparison** — after stripping the DP/PO prefixes, the
   pure numeric IDs are compared; the intersection is empty

3. **Demographic cross-comparison** — exact matching on the (age, sex, BMI)
   triplet. If the same individual appeared in both cohorts, the triplet would
   be identical; this is used to screen for potentially shared donors

---

## Note on sample counts

This check spans the **full GSE76895 series (103 samples)**: 32 non-diabetic,
36 T2D, 15 impaired-glucose-tolerance (IGT) and 20 type-3c-diabetes (T3cD)
donors. The manuscript's differential-expression analyses used only the
non-diabetic and T2D subset of GSE76895 (n = 68). The overlap check deliberately
covers all 103 donors, because a donor-overlap assessment must consider every
donor in each series, not only the subset retained for downstream analysis.

---

## Usage

```bash
python donor_overlap_check.py <GSE76895_series_matrix.txt.gz> <GSE76894_series_matrix.txt.gz>
```

With no arguments, the script reads the two series-matrix files under
`geo_data/` by default.

---

## Results

```
GSE76895 : 103 samples, 103 unique donor codes (DP, numbers 2–197)
GSE76894 : 103 samples,  82 unique donor codes (PO, numbers 755–946)

Layer 1  donor identifier intersection: 0
Layer 2  numeric-ID intersection: 0
Layer 3  exact (age,sex,BMI) triplet match: 0

Conclusion: no donor overlap between the two cohorts.
```

**Key facts**:
- GSE76895 consists of **partial pancreatectomy patients** (partial pancreatectomy; donor code DP);
- GSE76894 consists of **organ donors** (organ donor; donor code PO);
- The two are distinct populations with different donor numbering systems; **no donor appears in both cohorts**.

---

## Limitations (disclosed in the response letter)

The GEO series-matrix metadata **does not provide a formal unique subject
identifier (e.g., a study-level subject ID or a de-identified ID)**; this check
relies on the donor codes embedded in the sample titles and on demographic
features. Therefore:
- We can confirm **no overlap in donor numbers and no exact demographic match**;
- but we **cannot formally exclude a very low-probability undetected overlap**
  (e.g., the same donor using completely different codes in the two cohorts with
  inconsistent metadata records — an extremely unlikely scenario that also
  contradicts the "PPP patients vs organ donors" cohort design itself).

This caveat is stated in the response letter and in the Limitations section of
the manuscript.
