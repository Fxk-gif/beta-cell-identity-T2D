"""Supplementary Figure 3. STRING protein-protein interaction network.

Protein-protein interaction network of the beta-cell identity genes built from
STRING v12 (combined_score >= 0.40). Of the 15 common-panel genes, 14 form a
densely connected network (67 edges, mean combined score 0.680). Node fill
color shows the validation-set log2FC (GSE76894); edge width is proportional to
the STRING combined score; gold rings mark the five highest-degree hubs
(INS, SLC2A2, NEUROD1, ABCC8, GCK).

Inputs:
  - ml_output/ppi_string_result.xlsx  : STRING edge table (combined_score)
  - ml_output/identity_log2fc_both.csv: per-gene log2FC (validation cohort)
"""

import csv
import math
import os
import re
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

from _paths import find_workspace

WS = find_workspace()
ML = os.path.join(WS, "ml_output")
OUT = os.path.join(WS, "figures")
os.makedirs(OUT, exist_ok=True)


def load_ppi_edges(path):
    """Read STRING edges from an xlsx export via low-level zip parsing."""
    with zipfile.ZipFile(path) as z:
        shared = z.read("xl/sharedStrings.xml").decode("utf-8")
        ss = re.findall(r"<t[^>]*>(.*?)</t>", shared, re.S)
        sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8")
    rows = re.findall(r"<row[^>]*>(.*?)</row>", sheet, re.S)

    def cell_value(c):
        t = re.search(r't="(\w+)"', c)
        v = re.search(r"<v>(.*?)</v>", c)
        if not v:
            return ""
        val = v.group(1)
        if t and t.group(1) == "s":
            return ss[int(val)]
        return val

    edges = []
    for row in rows[1:]:  # skip header
        cells = re.findall(r"<c[^>]*>.*?</c>", row, re.S)
        vals = [cell_value(c) for c in cells]
        if len(vals) >= 12:
            edges.append((vals[0], vals[1], float(vals[11])))
    return edges


edges = load_ppi_edges(os.path.join(ML, "ppi_string_result.xlsx"))

# validation-set log2FC
fc = {}
with open(os.path.join(ML, "identity_log2fc_both.csv"),
          encoding="utf-8") as f:
    for r in csv.DictReader(f):
        fc[r["gene"]] = float(r["log2FC_GSE76894"])

# ---- circular layout ---------------------------------------------------------
nodes = sorted({e[0] for e in edges} | {e[1] for e in edges})
tf_group = ["HNF1A", "HNF4A", "NEUROD1", "PAX6", "MAFB"]
func_group = [g for g in nodes if g not in tf_group]
all_nodes = tf_group + func_group

pos = {}
for i, g in enumerate(all_nodes):
    ang = 2 * math.pi * i / len(all_nodes) - math.pi / 2
    pos[g] = (math.cos(ang), math.sin(ang))

# ---- figure ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 9), dpi=300)

for a, b, s in edges:
    x1, y1 = pos[a]
    x2, y2 = pos[b]
    lw = 0.5 + 3.0 * (s - 0.40) / (0.99 - 0.40)
    alpha = 0.35 + 0.65 * (s - 0.40) / (0.99 - 0.40)
    ax.plot([x1, x2], [y1, y2], color="#7f7f7f", lw=lw, alpha=alpha, zorder=1)

norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.1)


def node_color(g):
    if g not in fc:
        return "#bdbdbd"
    return plt.cm.RdBu_r(norm(fc[g]))


hub = {"INS", "SLC2A2", "NEUROD1", "ABCC8", "GCK"}
for g in all_nodes:
    x, y = pos[g]
    ax.add_patch(Circle((x, y), 0.115, fc=node_color(g),
                        ec="black", lw=1.2, zorder=3))
    if g in hub:
        ax.add_patch(Circle((x, y), 0.140, fc="none",
                            ec="#d4a017", lw=2.0, zorder=2))
    ax.text(x * 1.32, y * 1.32, g, ha="center", va="center", fontsize=9,
            fontweight="bold", color="#222222", zorder=4)

ax.set_xlim(-1.7, 1.7)
ax.set_ylim(-1.7, 1.7)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("STRING protein\u2013protein interaction network of the "
             "\u03b2-cell identity program\n"
             "14 identity genes form a densely connected network "
             "(67 edges, mean combined score = 0.680)",
             fontsize=12, fontweight="bold", pad=10)

handles = [
    Line2D([0], [0], color="#7f7f7f", lw=1,
           label="edge width \u221d STRING combined score"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#2b6cb0",
           markersize=10, label="downregulated (log\u2082FC < 0)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#c53030",
           markersize=10, label="upregulated (log\u2082FC > 0)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#bdbdbd",
           markersize=10, label="no probe"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
           markeredgecolor="#d4a017", markeredgewidth=2, markersize=12,
           label="network hub (top degree)"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True,
          framealpha=0.9, bbox_to_anchor=(0.0, 1.02))

fig.savefig(os.path.join(OUT, "SupplementaryFigure3.TIF"), dpi=300,
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

print("nodes:", len(nodes), "edges:", len(edges))
print("mean combined score:",
      round(sum(e[2] for e in edges) / len(edges), 3))
print("Supplementary Figure 3 written to", OUT)
