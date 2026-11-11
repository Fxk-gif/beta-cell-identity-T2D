"""Figure 6. Cross-dataset concordance of single-cell identity gene changes.

Panel (A): scatter of donor-level log2 fold-change (T2D vs ND) in the Xin
(GSE81608) and Baron (GSE84133) single-cell datasets, one point per identity
gene. Panel (B): the genes reduced in both datasets (SLC2A2, IAPP, INS),
shown as directional evidence only.

The donor-level statistics are read from the result tables produced by the
donor-level re-analysis scripts (xin_donor_level.py and baron_donor_level.py):

  - data/xin_donor_level_results.csv  : Xin donor-level Mann-Whitney U + BH
  - data/baron_donor_level_log2fc.csv : Baron donor-level log2FC (descriptive)

Xin has 18 donors (12 ND, 6 T2D), so a donor-level test is valid. Baron has a
single T2D donor, so its changes are directional only (no test is possible).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _paths import find_workspace

WS = find_workspace()
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(WS, "figures")
os.makedirs(OUT, exist_ok=True)

ND_C = "#4DBBD5"
T2D_C = "#E64B35"

# ---- load donor-level results ------------------------------------------------
xin = pd.read_csv(os.path.join(DATA, "xin_donor_level_results.csv"))
baron = pd.read_csv(os.path.join(DATA, "baron_donor_level_log2fc.csv"))

m = xin.merge(baron[["gene", "log2FC"]], on="gene",
              suffixes=("_xin", "_baron"))
m = m.dropna(subset=["log2FC_xin", "log2FC_baron"])

# ---- significance flag (Xin donor-level, BH) --------------------------------
m["sig"] = np.where((m["FDR"] < 0.05) & (m["log2FC_xin"] < 0), "down",
            np.where((m["FDR"] < 0.05) & (m["log2FC_xin"] > 0), "up", "ns"))

# ---- panel A: scatter -------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.5, 4.2), dpi=300)

color = {"down": T2D_C, "up": "#3B6FB6", "ns": "#9E9E9E"}
for _, r in m.iterrows():
    axA.scatter(r["log2FC_xin"], r["log2FC_baron"],
                color=color[r["sig"]], s=60, zorder=3)
    axA.annotate(r["gene"], (r["log2FC_xin"], r["log2FC_baron"]),
                 textcoords="offset points", xytext=(5, 4), fontsize=7)

axA.axhline(0, color="#b3b3b3", linewidth=0.8)
axA.axvline(0, color="#b3b3b3", linewidth=0.8)
axA.set_xlabel("log2FC (Xin, donor level)")
axA.set_ylabel("log2FC (Baron, donor level)")
axA.set_title("(A) Cross-dataset concordance", fontsize=10)
axA.set_aspect("equal", adjustable="box")

# ---- panel B: genes reduced in both -----------------------------------------
# The three genes reduced in both single-cell datasets and concordant with the
# bulk signal (SLC2A2, IAPP, INS), presented as directional evidence only.
CORE = ["SLC2A2", "IAPP", "INS"]
shared = m[m["gene"].isin(CORE)].copy()
shared = shared.sort_values("log2FC_xin")
bar_x = np.arange(len(shared))
width = 0.38
axB.barh(bar_x + width / 2, shared["log2FC_xin"], width,
         color=T2D_C, label="Xin")
axB.barh(bar_x - width / 2, shared["log2FC_baron"], width,
         color=ND_C, label="Baron")
axB.set_yticks(bar_x)
axB.set_yticklabels(shared["gene"])
axB.axvline(0, color="#808080", linewidth=0.8)
axB.set_xlabel("log2FC (T2D vs ND)")
axB.set_title("(B) Genes reduced in both datasets", fontsize=10)
axB.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure6.TIF"), dpi=300,
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

print("Genes reduced in both datasets:", ", ".join(shared["gene"].tolist()))
print("Figure 6 written to", OUT)
