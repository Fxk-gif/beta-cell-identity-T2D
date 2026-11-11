#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
overlap_corrected_gsea.py
=========================

Quantifies the overlap between the curated PDX1/HNF1A target gene sets and the
16-gene beta-cell identity score, and re-runs the GSEA after removing the
identity-score genes from each curated target set.

The script performs three steps:

  (1) Quantifies the overlap between the curated PDX1/HNF1A target gene sets
      and the 16-gene beta-cell identity score.
  (2) Re-runs the GSEA after removing all identity-score genes from each
      curated target set (de-circularized GSEA), using the weighted
      running-sum ES algorithm with seed = 42 and 2000 permutations.
  (3) Reports whether PDX1 appears in the HNF1A target list, and whether HNF1A
      appears in the PDX1 target list.

Inputs:
  - r3_inputs/volcano_input.csv     : whole-genome log2FC (discovery cohort)
  - r3_inputs/tf_target_sets.csv    : curated TF target sets

The 16-gene identity score (PDX1, NKX6-1, HNF1A, HNF4A, MAFB, NEUROD1, PAX6,
SLC2A2, GCK, ABCC8, IAPP, GLP1R, SLC30A8, PCSK1, PCSK2, INS).

Requires: numpy.

Usage:
    python overlap_corrected_gsea.py
"""

import csv
import os
import numpy as np

# ---- paths (resolve relative to this script's location) -----------------
HERE = os.path.dirname(os.path.abspath(__file__))
# The data reside in the analysis workspace r3_inputs/ directory. The location is
# resolved via the T2D_R3_INPUTS environment variable; if unset, the script walks
# upward to locate the project root containing r3_inputs/.
def _find_r3_inputs():
    env = os.environ.get("T2D_R3_INPUTS", "").rstrip("/\\")
    if env and os.path.isdir(env):
        return env
    d = os.path.abspath(HERE)
    for _ in range(12):
        if os.path.isdir(os.path.join(d, "r3_inputs")):
            return os.path.join(d, "r3_inputs")
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return os.path.join(HERE, "..", "..", "r3_inputs")

DATA_DIR = _find_r3_inputs()

VOLCANO = os.path.join(DATA_DIR, "volcano_input.csv")
TFSETS  = os.path.join(DATA_DIR, "tf_target_sets.csv")

# ---- the 16-gene identity score ------------------------------------------
IDENTITY_16 = {
    # master transcription factors (7)
    "PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
    # functional genes (9)
    "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8", "PCSK1", "PCSK2", "INS",
}

# ---- load data -----------------------------------------------------------
fc = {}
with open(VOLCANO, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        fc[row["id"]] = float(row["logFC"])

sets = {}
with open(TFSETS, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        sets.setdefault(row["set"], []).append(row["gene"])

# ---- GSEA (weighted running-sum ES) --------------------------------------
ranked = sorted(fc.items(), key=lambda x: -x[1])   # log2FC descending
N = len(ranked)


def compute_es(gene_set, ranked, fc, n_miss, nr):
    """Weighted running-sum ES. `nr` is the FIXED normalizer (the observed
    gene set's total |logFC|), held constant across observed and permutation
    samples (standard GSEA weighting)."""
    n_hits = len(gene_set)
    if nr == 0:
        return 0.0
    es = 0.0
    es_max = 0.0
    for g, lfc in ranked:
        if g in gene_set:
            es += abs(lfc) / nr
        else:
            es -= 1.0 / n_miss
        if abs(es) > abs(es_max):
            es_max = es
    return es_max


def gsea(geneset, n_perm=2000, seed=42):
    geneset = set(g for g in geneset if g in fc)
    if len(geneset) < 3:
        return None, None, None
    n_hits = len(geneset)
    n_miss = N - n_hits
    nr = sum(abs(fc[g]) for g in geneset)  # fixed normalizer from observed set
    es_max = compute_es(geneset, ranked, fc, n_miss, nr)
    all_genes = list(fc.keys())
    rng = np.random.default_rng(seed)
    perm_es = np.empty(n_perm)
    for i in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es[i] = compute_es(sample, ranked, fc, n_miss, nr)
    p = float(np.mean(np.abs(perm_es) >= abs(es_max)))
    return es_max, p, n_hits


# ---- (1) overlap quantification ------------------------------------------
print("=" * 80)
print("(1) OVERLAP between curated target sets and the 16-gene identity score")
print("=" * 80)
print(f"Identity-score genes (n={len(IDENTITY_16)}): "
      f"{', '.join(sorted(IDENTITY_16))}")
print()

for sname, genes in sets.items():
    gset = set(genes)
    shared = gset & IDENTITY_16
    print(f"  {sname:24s} size={len(gset):2d}  shared_with_score={len(shared):2d}  "
          f"[{', '.join(sorted(shared)) if shared else 'none'}]")

print()

# ---- (2) overlap-corrected GSEA ------------------------------------------
print("=" * 80)
print("(2) OVERLAP-CORRECTED GSEA (remove identity-score genes, re-run)")
print("=" * 80)

print("\n  -- sanity check: original GSEA results --")
print(f"  {'set':24s} {'ES':>9s} {'p':>8s} {'n_hits':>6s}")
orig = {}
for sname, genes in sets.items():
    es, p, n_hits = gsea(genes)
    orig[sname] = (es, p, n_hits)
    if es is not None:
        print(f"  {sname:24s} {es:+9.4f} {p:8.4f} {n_hits:6d}")

print("\n  -- overlap-corrected (de-circularized) GSEA --")
print(f"  {'set':24s} {'kept':>5s} {'ES':>9s} {'p':>8s}  verdict")
corrected = {}
for sname, genes in sets.items():
    gset = set(g for g in genes if g in fc)
    kept = [g for g in gset if g not in IDENTITY_16]
    es, p, n_hits = gsea(kept)
    corrected[sname] = (es, p, n_hits, len(kept))
    if es is None:
        print(f"  {sname:24s} {len(kept):5d}   {'--':>9s} {'--':>8s}  "
              f"<3 genes remain, cannot test")
        continue
    sig = "significant" if p < 0.05 else "NOT significant"
    direction = "down" if es < 0 else "up"
    print(f"  {sname:24s} {len(kept):5d} {es:+9.4f} {p:8.4f}  {sig} ({direction})")

# ---- (3) curation-fact verification --------------------------------------
print()
print("=" * 80)
print("(3) Curation-fact verification")
print("=" * 80)
pdx = set(sets.get("PDX1_targets", []))
hnf = set(sets.get("HNF1A_targets", []))
print(f"  PDX1 in HNF1A_targets list?  {'YES' if 'PDX1' in hnf else 'NO'}")
print(f"  HNF1A in PDX1_targets list?  {'YES' if 'HNF1A' in pdx else 'NO'}")
print(f"  HNF1A_targets (n={len(hnf)}): {', '.join(sorted(hnf))}")
print(f"  PDX1_targets  (n={len(pdx)}): {', '.join(sorted(pdx))}")
print()

# ---- summary table to CSV ------------------------------------------------
out = os.path.join(HERE, "overlap_corrected_gsea_results.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["set", "orig_ES", "orig_p", "orig_n_hits",
                "removed_identity_genes", "kept_n", "corrected_ES", "corrected_p"])
    for sname, genes in sets.items():
        gset = set(genes)
        shared = sorted(gset & IDENTITY_16)
        oe, op, onh = orig[sname]
        ce, cp, cn, ckept = corrected[sname]
        w.writerow([sname,
                    f"{oe:.6f}" if oe is not None else "",
                    f"{op:.6f}" if op is not None else "",
                    onh if onh is not None else "",
                    ";".join(shared),
                    ckept,
                    f"{ce:.6f}" if ce is not None else "",
                    f"{cp:.6f}" if cp is not None else ""])
print("=" * 80)
