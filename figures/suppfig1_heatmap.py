"""Supplementary Figure 1. Expression heatmap of beta-cell identity genes.

Group-mean expression (z-scored per gene) of the identity-related genes across
three categories -- beta-cell identity genes (master transcription factors and
functional genes), progenitor/disallowed genes, and neuroendocrine markers --
in the discovery (GSE76895) and validation (GSE76894) cohorts.

Category labels follow identity_gene_sets.csv, with the endocrine-identity and
neuroendocrine-transdifferentiation sets combined into a single neuroendocrine
category to match the figure legend.
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

# ---- category mapping --------------------------------------------------------
CATEGORY = {
    "beta_identity": "Beta-cell identity",
    "dediff": "Progenitor / disallowed",
    "endocrine_identity": "Neuroendocrine markers",
    "neuroendocrine_transdiff": "Neuroendocrine markers",
}
CAT_ORDER = ["Beta-cell identity", "Progenitor / disallowed",
             "Neuroendocrine markers"]

sets = pd.read_csv(os.path.join(DATA, "identity_gene_sets.csv"),
                   encoding="utf-8-sig")
gene_cat = {r["gene"]: CATEGORY[r["category"]] for _, r in sets.iterrows()}


def load(ds):
    expr = pd.read_csv(os.path.join(DATA, f"{ds}_expr.csv"), index_col=0)
    grp = pd.read_csv(os.path.join(DATA, f"{ds}_groups.csv"), header=None)
    grp.columns = ["sample", "group"]
    grp = grp[grp["group"].isin(["T2D", "ND"])]
    expr = expr[grp["sample"]]
    return expr, grp


def group_mean_z(ds):
    expr, grp = load(ds)
    genes = [g for g in gene_cat if g in expr.index]
    sub = expr.loc[genes].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    nd = z.loc[:, grp["group"].values == "ND"].mean(axis=1)
    t2d = z.loc[:, grp["group"].values == "T2D"].mean(axis=1)
    mat = pd.DataFrame({"ND": nd.values, "T2D": t2d.values},
                       index=genes)
    mat["category"] = [gene_cat[g] for g in mat.index]
    return mat


# ---- build combined matrix aligned across cohorts ---------------------------
m95 = group_mean_z("GSE76895")
m94 = group_mean_z("GSE76894")

# order genes by category then ND value (discovery cohort)
m95["cat_rank"] = m95["category"].map({c: i for i, c in enumerate(CAT_ORDER)})
m95 = m95.sort_values(["cat_rank", "ND"])
ordered_genes = m95.index.tolist()
# keep the same gene order for the validation cohort
m94 = m94.loc[[g for g in ordered_genes if g in m94.index]]

# ---- figure: one heatmap per cohort, stacked vertically ---------------------
fig, axes = plt.subplots(2, 1, figsize=(3.6, 9.5), dpi=300,
                         gridspec_kw={"height_ratios":
                                      [len(m95), len(m94)]})

for ax, (mat, ds) in zip(axes, [(m95, "GSE76895"), (m94, "GSE76894")]):
    data = mat[["ND", "T2D"]].to_numpy(dtype=float)
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                   interpolation="nearest")
    ax.set_yticks(np.arange(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["ND", "T2D"], fontsize=8)
    ax.set_title(ds, fontsize=10, pad=6)
    # category separators (drawn on the left as grouped labels)
    prev = None
    for i, c in enumerate(mat["category"]):
        if c != prev:
            if prev is not None:
                ax.axhline(i - 0.5, color="black", linewidth=0.8)
            prev = c
    # tick color by category
    for tick, c in zip(ax.get_yticklabels(), mat["category"]):
        tick.set_color("#1a1a1a")
    ax.grid(False)
    # draw category labels on the right side
    for c in CAT_ORDER:
        idx = mat.index[mat["category"] == c]
        if len(idx):
            pos = mat.index.get_loc(idx[-1]) + 0.5
            ax.text(2.25, pos - 0.5, c, va="center", ha="left",
                    fontsize=7, fontweight="bold", rotation=0)

# shared colorbar (vertical, on the right)
fig.subplots_adjust(left=0.22, right=0.86, top=0.90, bottom=0.03, hspace=0.35)
cbar_ax = fig.add_axes([0.90, 0.25, 0.02, 0.5])
fig.colorbar(im, cax=cbar_ax, label="Group-mean z-scored expression")

fig.suptitle("Group-mean expression of identity-related genes",
             fontsize=11, y=0.99)
fig.savefig(os.path.join(OUT, "SupplementaryFigure1.TIF"), dpi=300,
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

print("Supplementary Figure 1 written to", OUT)
