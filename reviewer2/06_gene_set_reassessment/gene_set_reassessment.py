# -*- coding: utf-8 -*-
"""
Response to Reviewer 2, Comment 6:
Reassessment of curated gene-set definitions, effective gene counts, and
enrichment statistics.

This script addresses three specific issues raised by the reviewer:
  1. Redundant identifiers (GLUT2/SLC2A2, PC1/PCSK1) are removed so that each
     gene set is defined strictly by official HGNC symbol with each gene
     appearing exactly once.
  2. The Figure 9 shared-target error is corrected: the intersection of the
     de-duplicated PDX1 and HNF1A target sets is three genes (SLC2A2, INS, GCK),
     not SLC2A2 alone.
  3. Genes lacking GPL570 (HG-U133 Plus 2.0) probe coverage (MAFA; GLUT2 as an
     SLC2A2 alias) are flagged and excluded from the effective gene count.

It then recomputes the effective gene count and GSEA enrichment statistics
(enrichment score and nominal p value) on the de-duplicated, coverage-filtered
sets for both the discovery (GSE76895) and validation (GSE76894) cohorts.

Important methodological note: the expression matrices are RMA-normalized and
already in log2 space, so the correct fold change is the difference of group
means in log2 space (m_T2D - m_ND), not log2(ratio). The discovery-cohort
volcano input already uses this convention; the validation-cohort log2FC is
recomputed from the expression matrix directly.
"""

import os
import csv
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Locate the data root (the directory that contains r3_inputs).
# ---------------------------------------------------------------------------
_here = os.path.abspath(os.path.dirname(__file__))
BASE = None
d = _here
for _ in range(10):
    if os.path.isdir(os.path.join(d, "r3_inputs")):
        BASE = d
        break
    parent = os.path.dirname(d)
    if parent == d:
        break
    d = parent

if BASE is None:
    BASE = os.environ.get("T2D_DATA", "")
    if not BASE or not os.path.isdir(os.path.join(BASE, "r3_inputs")):
        raise FileNotFoundError(
            "Could not locate the data root containing r3_inputs. "
            "Set the T2D_DATA environment variable to the project directory."
        )

R3IN = os.path.join(BASE, "r3_inputs")
GEO = os.path.join(BASE, "geo_data", "processed")

# ---------------------------------------------------------------------------
# Gene-set definitions before de-duplication (these lists contain the
# redundant legacy-name entries that are removed below).
# ---------------------------------------------------------------------------
# Legacy-name alias map: old name -> official HGNC symbol.
ALIAS = {
    "GLUT2": "SLC2A2",   # GLUT2 is the historical name of SLC2A2
    "PC1": "PCSK1",      # PC1 is the legacy name of PCSK1
}

# Original PDX1 target set (contains the GLUT2/SLC2A2 redundancy).
PDX1_TARGETS_RAW = [
    "INS", "IAPP", "GCK", "SLC2A2", "NKX6-1", "MAFA", "PAX4",
    "GLUT2", "PC", "G6PC2",
]

# Original HNF1A target set (no redundancy in this set).
HNF1A_TARGETS_RAW = [
    "SLC2A2", "PKLR", "HNF4A", "HNF4G", "FOXA3", "INS", "GCK",
    "SLC30A8", "PDX1", "ABCC8", "KCNJ11", "HNF1B",
]

# Original beta-cell functional gene set (contains the PC1/PCSK1 redundancy).
BETA_FUNCTIONAL_RAW = [
    "INS", "IAPP", "SLC2A2", "GCK", "ABCC8", "KCNJ11", "SLC30A8",
    "PC", "PC1", "PCSK1", "PCSK2", "G6PC2", "GLP1R", "UCN3", "CHGA", "SYP",
]

# ---------------------------------------------------------------------------
# De-duplication: map legacy names to official symbols and collapse duplicates.
# ---------------------------------------------------------------------------
def deduplicate(genes):
    """Return the de-duplicated list of official gene symbols, preserving order."""
    seen = []
    for g in genes:
        official = ALIAS.get(g, g)
        if official not in seen:
            seen.append(official)
    return seen


PDX1_TARGETS = deduplicate(PDX1_TARGETS_RAW)
HNF1A_TARGETS = deduplicate(HNF1A_TARGETS_RAW)
BETA_FUNCTIONAL = deduplicate(BETA_FUNCTIONAL_RAW)


# ---------------------------------------------------------------------------
# GPL570 probe-coverage map.
# ---------------------------------------------------------------------------
def load_probe_coverage():
    sym2probes = {}
    path = os.path.join(GEO, "gpl570_probe2symbol.csv")
    with open(path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        next(rdr)  # header
        for r in rdr:
            probe, sym = r[0], r[1]
            sym2probes.setdefault(sym, []).append(probe)
    return sym2probes


def coverage_table(sym2probes, genes):
    """Return dict gene -> number of GPL570 probes (0 means no coverage)."""
    return {g: len(sym2probes.get(g, [])) for g in genes}


# ---------------------------------------------------------------------------
# Expression loading and corrected log2FC.
# ---------------------------------------------------------------------------
def load_expr(ds):
    expr_path = os.path.join(R3IN, f"{ds}_expr.csv")
    grp_path = os.path.join(R3IN, f"{ds}_groups.csv")
    with open(expr_path, encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = [r for r in rdr]
    genes = [r[0] for r in rows]
    mat = np.array([[float(x) for x in r[1:]] for r in rows])
    gmap = {}
    with open(grp_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gmap[row["sample"]] = row["group"]
    return genes, header[1:], mat, gmap


def full_log2fc(ds):
    """Corrected per-gene log2FC for ALL genes (m_T2D - m_ND, log2 space)."""
    genes, samples, mat, gmap = load_expr(ds)
    t2d_idx = [j for j, s in enumerate(samples) if gmap.get(s) == "T2D"]
    nd_idx = [j for j, s in enumerate(samples) if gmap.get(s) == "ND"]
    fc = {}
    for gi, g in enumerate(genes):
        m_t = mat[gi, t2d_idx].mean()
        m_n = mat[gi, nd_idx].mean()
        fc[g] = m_t - m_n
    return fc


def load_discovery_logfc():
    """Discovery-cohort (GSE76895) whole-genome logFC from volcano_input.csv."""
    fc = {}
    path = os.path.join(R3IN, "volcano_input.csv")
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            fc[row["id"]] = float(row["logFC"])
    return fc


# ---------------------------------------------------------------------------
# GSEA (single-tailed enrichment of a gene set within a ranked list).
# ---------------------------------------------------------------------------
def gsea(geneset, ranked, fc, n_perm=2000, seed=42):
    """
    Compute the enrichment score (ES) and nominal p value for a gene set
    within a gene list ranked by descending logFC, using a permutation test.
    geneset: list of official symbols with >=1 probe coverage.
    ranked: list of (gene, logFC) sorted descending by logFC.
    """
    gset = [g for g in geneset if g in fc]
    if len(gset) < 3:
        return None, None, None, None
    n_hits = len(gset)
    N = len(ranked)
    n_miss = N - n_hits
    # weights: |logFC| of hits
    nr = sum(abs(fc[g]) for g in gset)
    if nr == 0:
        return None, None, None, None

    def es_of(hit_set):
        es = 0.0
        es_max = 0.0
        for g, lfc in ranked:
            if g in hit_set:
                es += abs(fc[g]) / nr
            else:
                es -= 1.0 / n_miss
            if abs(es) > abs(es_max):
                es_max = es
        return es_max

    es_obs = es_of(set(gset))
    all_genes = list(fc.keys())
    rng = np.random.default_rng(seed)
    perm_es = []
    for _ in range(n_perm):
        sample = set(rng.choice(all_genes, size=n_hits, replace=False))
        perm_es.append(es_of(sample))
    perm_es = np.array(perm_es)
    p_nom = float(np.mean(np.abs(perm_es) >= abs(es_obs)))
    return es_obs, p_nom, n_hits, n_miss


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def main():
    sym2probes = load_probe_coverage()

    print("=" * 78)
    print("1. De-duplication of gene sets (legacy-name alias removal)")
    print("=" * 78)
    print(f"PDX1 targets:    raw {len(PDX1_TARGETS_RAW)} -> "
          f"deduplicated {len(PDX1_TARGETS)}")
    print(f"  removed: {[g for g in PDX1_TARGETS_RAW if g in ALIAS]}")
    print(f"HNF1A targets:   raw {len(HNF1A_TARGETS_RAW)} -> "
          f"deduplicated {len(HNF1A_TARGETS)}")
    print(f"beta functional: raw {len(BETA_FUNCTIONAL_RAW)} -> "
          f"deduplicated {len(BETA_FUNCTIONAL)}")
    print(f"  removed: {[g for g in BETA_FUNCTIONAL_RAW if g in ALIAS]}")

    print()
    print("=" * 78)
    print("2. Shared PDX1/HNF1A targets (corrected Figure 9)")
    print("=" * 78)
    shared = sorted(set(PDX1_TARGETS) & set(HNF1A_TARGETS))
    print(f"PDX1 targets:    {sorted(PDX1_TARGETS)}")
    print(f"HNF1A targets:   {sorted(HNF1A_TARGETS)}")
    print(f"Intersection:    {shared}")
    print(f"Number shared:   {len(shared)} (not SLC2A2 alone)")

    print()
    print("=" * 78)
    print("3. GPL570 probe coverage of every gene in the de-duplicated sets")
    print("=" * 78)
    all_genes = sorted(set(PDX1_TARGETS + HNF1A_TARGETS + BETA_FUNCTIONAL))
    cov = coverage_table(sym2probes, all_genes)
    no_cov = []
    print(f"{'gene':10s} {'probes':>7s}  {'set(s)':s}")
    set_map = {
        g: ("P" if g in PDX1_TARGETS else "") +
           ("H" if g in HNF1A_TARGETS else "") +
           ("F" if g in BETA_FUNCTIONAL else "")
        for g in all_genes
    }
    for g in all_genes:
        n = cov[g]
        tag = "  <-- NO COVERAGE" if n == 0 else ""
        if n == 0:
            no_cov.append(g)
        print(f"{g:10s} {n:7d}  {set_map[g]}{tag}")
    print(f"\nGenes without GPL570 coverage: {no_cov}")

    print()
    print("=" * 78)
    print("4. Effective gene counts (deduplicated + coverage-filtered)")
    print("=" * 78)
    sets = {
        "PDX1_targets": PDX1_TARGETS,
        "HNF1A_targets": HNF1A_TARGETS,
        "beta_functional": BETA_FUNCTIONAL,
    }
    eff_count = {}
    for sname, genes in sets.items():
        listed = len(genes)
        covered = [g for g in genes if cov[g] > 0]
        eff_count[sname] = covered
        print(f"{sname:18s}: listed {listed:2d} -> effective {len(covered):2d} "
              f"(excluded: {[g for g in genes if cov[g] == 0]})")

    print()
    print("=" * 78)
    print("5. GSEA enrichment statistics (corrected log2FC, both cohorts)")
    print("=" * 78)

    # Discovery cohort: whole-genome logFC from volcano_input.csv (correct).
    fc_disc = load_discovery_logfc()
    # Validation cohort: recompute whole-genome corrected log2FC from matrix.
    fc_val = full_log2fc("GSE76894")

    for cohort_name, fc in [("GSE76895 (discovery)", fc_disc),
                            ("GSE76894 (validation)", fc_val)]:
        ranked = sorted(fc.items(), key=lambda x: -x[1])
        print(f"\n--- {cohort_name} ---")
        print(f"  Whole-genome genes with logFC: {len(fc)}")
        print(f"{'gene set':18s} {'eff_n':>6s} {'ES':>8s} {'p_nominal':>10s}")
        results = {}
        for sname in ["PDX1_targets", "HNF1A_targets", "beta_functional"]:
            covered = eff_count[sname]
            es, p, n_hits, n_miss = gsea(covered, ranked, fc)
            if es is None:
                print(f"{sname:18s} {'--':>6s} {'--':>8s} {'--':>10s}  (too few genes)")
                continue
            results[sname] = (es, p, n_hits)
            print(f"{sname:18s} {n_hits:6d} {es:+8.4f} {p:10.4f}")

    # Per-gene membership, probe coverage, and effective-count table.
    out_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "gene_set_reassessment.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "PDX1_target", "HNF1A_target", "beta_functional",
                    "GPL570_probes", "covered"])
        for g in all_genes:
            w.writerow([
                g,
                "yes" if g in PDX1_TARGETS else "",
                "yes" if g in HNF1A_TARGETS else "",
                "yes" if g in BETA_FUNCTIONAL else "",
                cov[g],
                "yes" if cov[g] > 0 else "no",
            ])

    print()
    print(f"Per-gene reassessment table: {os.path.basename(out_csv)}")


if __name__ == "__main__":
    main()
