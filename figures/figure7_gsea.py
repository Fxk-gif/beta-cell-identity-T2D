"""Figure 7. Gene set enrichment analysis of transcription-factor target genes.

Panel (A): GSEA running-enrichment curves for the four curated target gene
sets in the original (unadjusted) analysis. Panel (B): the same curves after
overlap correction (removing genes shared with the 16-gene identity score),
which shows that only the beta-cell functional gene set retains significant
enrichment after de-circularization.

The weighted running-sum enrichment score is computed exactly as in
overlap_corrected_gsea.py (seed = 42, 2000 permutations).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _paths import find_workspace

WS = find_workspace()
DATA = os.path.join(WS, "r3_inputs")
OUT = os.path.join(WS, "figures")
os.makedirs(OUT, exist_ok=True)

# ---- 16-gene identity score (same definition as the GSEA correction script) --
IDENTITY_16 = {
    "PDX1", "NKX6-1", "HNF1A", "HNF4A", "MAFB", "NEUROD1", "PAX6",
    "SLC2A2", "GCK", "ABCC8", "IAPP", "GLP1R", "SLC30A8", "PCSK1", "PCSK2", "INS",
}

# ---- load ranked gene list and target sets -----------------------------------
fc = pd.read_csv(os.path.join(DATA, "volcano_input.csv"), encoding="utf-8-sig")
fc = fc[["id", "logFC"]].dropna()
ranked = fc.sort_values("logFC", ascending=False).reset_index(drop=True)
ranked_list = list(zip(ranked["id"], ranked["logFC"]))

sets = {}
tf = pd.read_csv(os.path.join(DATA, "tf_target_sets.csv"), encoding="utf-8-sig")
for _, r in tf.iterrows():
    sets.setdefault(r["set"], []).append(r["gene"])

N = len(ranked_list)
fc_map = dict(ranked_list)


def running_es(geneset, ranked_list, nr):
    """Weighted running-sum enrichment score (observed or permutation)."""
    n_hits = len(geneset)
    n_miss = N - n_hits
    if nr == 0 or n_miss == 0:
        return np.zeros(N)
    es = np.zeros(N)
    running = 0.0
    for i, (g, lfc) in enumerate(ranked_list):
        if g in geneset:
            running += abs(lfc) / nr
        else:
            running -= 1.0 / n_miss
        es[i] = running
    return es


def gsea_curve(geneset):
    geneset = [g for g in geneset if g in fc_map]
    if len(geneset) < 3:
        return None, None, None
    nr = sum(abs(fc_map[g]) for g in geneset)
    es = running_es(set(geneset), ranked_list, nr)
    es_max = max(es.min(), es.max(), key=abs)
    return es, es_max, len(geneset)


# ---- compute curves for all four sets (original and corrected) ---------------
SET_ORDER = ["PDX1_targets", "HNF1A_targets",
             "beta_functional_genes", "beta_identity_network"]
SET_LABEL = {
    "PDX1_targets": "PDX1 targets",
    "HNF1A_targets": "HNF1A targets",
    "beta_functional_genes": "beta-cell functional genes",
    "beta_identity_network": "beta-cell identity network",
}
COLORS = ["#E64B35", "#3B6FB6", "#2F9E6F", "#C792EA"]

fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), dpi=300)
fig.suptitle("GSEA: TF target gene sets enriched among downregulated genes",
             fontsize=12, y=0.98)

for ax, sname, color in zip(axes.flat, SET_ORDER, COLORS):
    genes = sets[sname]

    # original (unadjusted)
    es_orig, es_max_orig, n_orig = gsea_curve(genes)
    # overlap-corrected
    kept = [g for g in genes if g not in IDENTITY_16 and g in fc_map]
    es_corr, es_max_corr, n_corr = gsea_curve(kept)

    if es_orig is not None:
        ax.plot(es_orig, color=color, linewidth=1.3,
                label=f"original (ES = {es_max_orig:+.3f})")
    if es_corr is not None:
        ax.plot(es_corr, color="#808080", linewidth=1.1, linestyle="--",
                label=f"overlap-corrected (ES = {es_max_corr:+.3f})")

    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_title(SET_LABEL[sname], fontsize=10)
    ax.set_xlabel("Rank (high to low logFC)")
    ax.set_ylabel("Enrichment score")
    ax.legend(fontsize=7, frameon=False, loc="lower left")

fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(OUT, "Figure7.TIF"), dpi=300,
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

# ---- report final numbers ----------------------------------------------------
print("GSEA enrichment scores (original / overlap-corrected):")
for sname in SET_ORDER:
    genes = sets[sname]
    es_o, mx_o, n_o = gsea_curve(genes)
    kept = [g for g in genes if g not in IDENTITY_16 and g in fc_map]
    es_c, mx_c, n_c = gsea_curve(kept)
    orig_str = f"orig ES={mx_o:+.3f} (n={n_o})" if mx_o is not None else "orig (n<3)"
    corr_str = (f"corrected ES={mx_c:+.3f} (n={n_c})"
                if mx_c is not None else "corrected (n<3)")
    print(f"  {sname:24s}  {orig_str}  {corr_str}")

print("Figure 7 written to", OUT)
