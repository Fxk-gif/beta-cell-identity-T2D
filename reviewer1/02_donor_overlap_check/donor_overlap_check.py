#!/usr/bin/env python3
"""
donor_overlap_check.py
=======================

Check for donor overlap between the two bulk microarray cohorts used in the
manuscript (Manuscript ID 2015413, Frontiers in Endocrinology):

  - GSE76895 : human pancreatic islets, partial-pancreatectomy (PPP) cohort,
               Affymetrix GPL570 (HG-U133 Plus 2.0)
  - GSE76894 : human pancreatic islets, organ-donor (OD) cohort,
               Affymetrix GPL570 (HG-U133 Plus 2.0)

Both datasets originate from the IMIDIA biobank and share the same platform,
which is why Reviewer 1 (Comment #2) raised the concern that the two cohorts
may not be truly "independent" and that donor overlap should be checked and
stated explicitly.

This script addresses that comment by performing a multi-layer donor-overlap
check on the publicly available GEO metadata:

  Layer 1  Donor identifier cross-match
           (sample titles carry donor IDs: "DP###" for PPP, "PO###" for OD)
  Layer 2  Numeric donor-code cross-match
           (strip the DP/PO prefix, compare the numeric codes)
  Layer 3  Demographic cross-match
           (age + sex + BMI triple match would flag a possible shared donor)

Note on sample counts: this check covers the full GSE76895 series (103 samples,
comprising 32 non-diabetic, 36 T2D, 15 IGT and 20 T3cD donors). The manuscript's
differential-expression analyses used only the non-diabetic and T2D subset of
GSE76895 (n = 68); the overlap check deliberately spans all 103 donors because a
donor-overlap assessment must consider every donor in each series, not only the
subset retained for downstream analysis.

Requirements: Python 3 (standard library only; no third-party packages).

Usage:
    python donor_overlap_check.py

The script reads the two GEO series-matrix files (gzipped) directly. Set the
paths below or pass them as command-line arguments:

    python donor_overlap_check.py <GSE76895_series_matrix.txt.gz> <GSE76894_series_matrix.txt.gz>
"""

import gzip
import os
import re
import sys
from collections import defaultdict

# --- default data paths -------------------------------------------------
# Locate the data root (the directory containing geo_data/) by walking up from
# this script, so the defaults remain valid when the script is relocated.
_here = os.path.abspath(os.path.dirname(__file__))
_data_root = None
_d = _here
for _ in range(10):
    if os.path.isdir(os.path.join(_d, "geo_data")):
        _data_root = _d
        break
    _parent = os.path.dirname(_d)
    if _parent == _d:
        break
    _d = _parent

if _data_root is None:
    _data_root = _here  # fall back to the script directory

DEFAULT_GSE76895 = os.path.join(_data_root, "geo_data", "GSE76895_series_matrix.txt.gz")
DEFAULT_GSE76894 = os.path.join(_data_root, "geo_data", "GSE76894_series_matrix.txt.gz")


def parse_series_matrix(path):
    """
    Parse a GEO SOFT series-matrix file (gzipped) and return, per sample:
        geo_id        -> GSM accession
        title         -> sample title (contains the donor code)
        characteristics -> dict of the sample characteristics
                          (age, sex, BMI, diabetes status, ...)
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        lines = fh.read().split("\n")

    def get_field(prefix):
        for line in lines:
            if line.startswith(prefix):
                return line.split("\t")[1:]
        return []

    geo_ids = get_field("!Sample_geo_accession")
    titles = get_field("!Sample_title")
    # characteristics appear on multiple repeated !Sample_characteristics_ch1 lines
    char_lines = [l.split("\t")[1:] for l in lines
                  if l.startswith("!Sample_characteristics_ch1")]

    n = len(geo_ids)
    samples = []
    for i in range(n):
        char = {}
        for cl in char_lines:
            if i < len(cl) and cl[i]:
                raw = cl[i].strip('"')
                if ":" in raw:
                    k, v = raw.split(":", 1)
                    char[k.strip()] = v.strip()
                else:
                    char[raw] = raw
        samples.append({
            "geo": geo_ids[i].strip('"') if i < len(geo_ids) else "",
            "title": titles[i].strip('"') if i < len(titles) else "",
            "char": char,
        })
    return samples


def extract_donor_code(title):
    """
    Extract the donor code from a GEO sample title.

    GSE76895 (PPP) titles look like:  H1_DP002PR01: Pancreatic islet ...
    GSE76894 (OD)  titles look like:  H1_PO755IR01: Pancreatic islet ...

    Returns (donor_prefix, donor_number) or (None, None) if not matched.
    """
    m = re.match(r"\"?(?:H\d+_)?(D[P]|PO)(\d+)", title)
    if m:
        return m.group(1), m.group(2)
    # fallback: any DPxxx / POxxx token anywhere in the title
    m2 = re.search(r"\b(DP|PO)(\d+)\b", title)
    if m2:
        return m2.group(1), m2.group(2)
    return None, None


def num(s):
    """Safe float conversion, returns None on failure."""
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    # resolve paths
    if len(sys.argv) >= 3:
        p95, p94 = sys.argv[1], sys.argv[2]
    else:
        p95, p94 = DEFAULT_GSE76895, DEFAULT_GSE76894

    print("=" * 78)
    print("Donor-overlap check: GSE76895 (PPP) vs GSE76894 (OD)")
    print("=" * 78)

    s95 = parse_series_matrix(p95)
    s94 = parse_series_matrix(p94)

    print(f"\nGSE76895 : {len(s95)} samples")
    print(f"GSE76894 : {len(s94)} samples")

    # ---- Layer 1 & 2 : donor identifier cross-match ---------------------
    donors95 = {}   # donor_code -> list of sample geos
    donors94 = {}
    for s in s95:
        pre, no = extract_donor_code(s["title"])
        if pre:
            donors95.setdefault(f"{pre}{no}", []).append(s["geo"])
    for s in s94:
        pre, no = extract_donor_code(s["title"])
        if pre:
            donors94.setdefault(f"{pre}{no}", []).append(s["geo"])

    print(f"\n[Layer 1] Donor identifier cross-match")
    print(f"  GSE76895 unique donor codes : {len(donors95)}  (prefix DP = pancreatectomy)")
    print(f"  GSE76894 unique donor codes : {len(donors94)}  (prefix PO = organ donor)")
    print(f"  Example GSE76895 donor codes: {sorted(donors95)[:6]}")
    print(f"  Example GSE76894 donor codes: {sorted(donors94)[:6]}")

    overlap_ids = set(donors95) & set(donors94)
    print(f"  Donor-code intersection      : {sorted(overlap_ids) if overlap_ids else 'NONE'}")
    print(f"  Overlapping donor codes      : {len(overlap_ids)}")

    print(f"\n[Layer 2] Numeric donor-code cross-match (strip DP/PO prefix)")
    nums95 = set(re.sub(r"\D", "", c) for c in donors95)
    nums94 = set(re.sub(r"\D", "", c) for c in donors94)
    overlap_nums = nums95 & nums94
    print(f"  GSE76895 numeric range       : {min(int(x) for x in nums95)} - {max(int(x) for x in nums95)}")
    print(f"  GSE76894 numeric range       : {min(int(x) for x in nums94)} - {max(int(x) for x in nums94)}")
    print(f"  Numeric-code intersection    : {sorted(int(x) for x in overlap_nums) if overlap_nums else 'NONE'}")
    print(f"  Overlapping numeric codes    : {len(overlap_nums)}")

    # ---- Layer 3 : demographic cross-match ------------------------------
    # Build (age, sex, BMI) triple keys; a triple match would flag a possible
    # shared donor (a real person cannot be in both PPP and OD cohorts).
    def demo_key(s):
        c = s["char"]
        age = num(re.sub(r"[^0-9.]", "", c.get("age (years)", "")))
        sex = c.get("gender (m/f)", "").strip().upper()
        bmi = num(re.sub(r"[^0-9.]", "", c.get("body mass index (bmi)", "")))
        # round BMI to 2dp and age to integer to tolerate minor formatting diffs
        return (None if age is None else round(age),
                sex or None,
                None if bmi is None else round(bmi, 2))

    demo95 = defaultdict(list)
    demo94 = defaultdict(list)
    for s in s95:
        demo95[demo_key(s)].append(s["geo"])
    for s in s94:
        demo94[demo_key(s)].append(s["geo"])

    print(f"\n[Layer 3] Demographic cross-match (age + sex + BMI exact triple)")
    # only consider keys where all three fields are present (non-None)
    usable95 = {k: v for k, v in demo95.items() if None not in k}
    usable94 = {k: v for k, v in demo94.items() if None not in k}
    print(f"  GSE76895 samples with complete (age,sex,BMI): {sum(len(v) for v in usable95.values())}")
    print(f"  GSE76894 samples with complete (age,sex,BMI): {sum(len(v) for v in usable94.values())}")
    shared_demo = set(usable95) & set(usable94)
    print(f"  Exact (age,sex,BMI) triples shared         : {len(shared_demo)}")
    if shared_demo:
        for k in sorted(shared_demo):
            print(f"    {k}: GSE76895={usable95[k]}  GSE76894={usable94[k]}")
    else:
        print("    (none)")

    # ---- conclusion ------------------------------------------------------
    print("\n" + "=" * 78)
    print("CONCLUSION")
    print("=" * 78)
    n_overlap = len(overlap_ids) + len(overlap_nums) + len(shared_demo)
    if n_overlap == 0:
        print("No donor overlap detected between GSE76895 (PPP) and GSE76894 (OD).")
        print("The two cohorts carry disjoint donor identifiers (DP### vs PO###),")
        print("disjoint numeric donor codes, and no exact (age, sex, BMI) triple")
        print("matches, consistent with two non-overlapping donor groups drawn from")
        print("distinct sample-collection protocols (partial pancreatectomy vs.")
        print("organ donation).")
    else:
        print("Potential overlap signal detected; review the flagged matches above.")
    print("=" * 78)


if __name__ == "__main__":
    main()
