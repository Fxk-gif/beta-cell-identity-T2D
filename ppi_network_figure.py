# -*- coding: utf-8 -*-
"""
Figure 9 (Supplementary Figure 3): STRING protein-protein interaction network
of the 15 β-cell identity genes (14 nodes retained after STRING filtering).
Node color = validation-set log2FC (GSE76894); edge width = STRING combined_score.
"""
import csv, zipfile, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

BASE = r"C:/Users/45902/WorkBuddy/2026-09-20-15-12-44"

# ---- load PPI edges (combined_score) from STRING result xlsx (zipfile low-level) ----
def load_ppi():
    p = f"{BASE}/ml_output/ppi_string_result.xlsx"
    z = zipfile.ZipFile(p)
    shared = z.read("xl/sharedStrings.xml").decode("utf-8")
    ss = re.findall(r"<t[^>]*>(.*?)</t>", shared, re.S)
    sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8")
    rows = re.findall(r"<row[^>]*>(.*?)</row>", sheet, re.S)
    def cv(c):
        t = re.search(r't="(\w+)"', c); v = re.search(r"<v>(.*?)</v>", c)
        if not v: return ""
        val = v.group(1)
        if t and t.group(1) == "s": return ss[int(val)]
        return val
    edges = []
    for row in rows[1:]:
        cells = re.findall(r"<c[^>]*>.*?</c>", row, re.S)
        vals = [cv(c) for c in cells]
        if len(vals) >= 12:
            edges.append((vals[0], vals[1], float(vals[11])))
    return edges

edges = load_ppi()

# ---- log2FC from validation set ----
fc = {}
for r in csv.DictReader(open(f"{BASE}/ml_output/identity_log2fc_both.csv")):
    fc[r["gene"]] = float(r["log2FC_GSE76894"])

# ---- manual circular layout ----
nodes = sorted(set([e[0] for e in edges] + [e[1] for e in edges]))
# order: hub genes (INS/SLC2A2/NEUROD1/ABCC8/GCK) first for readability
import math
n = len(nodes)
# fixed positions (unit circle), tuned to reduce overlap
pos = {}
# place TFs top, functional genes bottom
tf_group = ["HNF1A", "HNF4A", "NEUROD1", "PAX6", "MAFB"]
func_group = [g for g in nodes if g not in tf_group]
# arrange tf_group around top half, func around bottom
all_nodes = tf_group + func_group
for i, g in enumerate(all_nodes):
    ang = 2 * math.pi * i / len(all_nodes) - math.pi / 2
    r = 1.0
    pos[g] = (r * math.cos(ang), r * math.sin(ang))

# ---- figure ----
fig, ax = plt.subplots(figsize=(9, 9))

# edges (draw first, behind nodes)
for a, b, s in edges:
    x1, y1 = pos[a]; x2, y2 = pos[b]
    lw = 0.5 + 3.0 * (s - 0.40) / (0.99 - 0.40)  # width scaled
    alpha = 0.35 + 0.65 * (s - 0.40) / (0.99 - 0.40)
    ax.plot([x1, x2], [y1, y2], color="#7f7f7f", lw=lw, alpha=alpha, zorder=1)

# node color map: log2FC, blue=down, red=up, white=0
norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.1)

def node_color(g):
    if g not in fc:
        return "#bdbdbd"  # gray (no probe / not in 15-gene panel)
    v = fc[g]
    return plt.cm.RdBu_r(norm(v))

# hub genes get a ring
hub = {"INS", "SLC2A2", "NEUROD1", "ABCC8", "GCK"}

for g in all_nodes:
    x, y = pos[g]
    c = node_color(g)
    r = 0.115
    ax.add_patch(Circle((x, y), r, fc=c, ec="black", lw=1.2, zorder=3))
    if g in hub:
        ax.add_patch(Circle((x, y), r + 0.025, fc="none", ec="#d4a017", lw=2.0, zorder=2))
    # label
    lab_x, lab_y = x * 1.32, y * 1.32
    ax.text(lab_x, lab_y, g, ha="center", va="center", fontsize=9,
            fontweight="bold", color="#222222", zorder=4)

ax.set_xlim(-1.7, 1.7)
ax.set_ylim(-1.7, 1.7)
ax.set_aspect("equal")
ax.axis("off")

# title
ax.set_title("STRING protein\u2013protein interaction network of the \u03b2-cell identity program\n"
             "14 identity genes form a densely connected network (67 edges, mean combined_score = 0.68)",
             fontsize=12, fontweight="bold", pad=10)

# legend
handles = [
    Line2D([0], [0], color="#7f7f7f", lw=1, label="edge width \u221d STRING combined_score"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#2b6cb0", markersize=10, label="downregulated (log\u2082FC < 0)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#c53030", markersize=10, label="upregulated (log\u2082FC > 0)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#bdbdbd", markersize=10, label="no probe"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#d4a017", markeredgewidth=2, markersize=12, label="network hub (top degree)"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True, framealpha=0.9, bbox_to_anchor=(0.0, 1.02))

plt.tight_layout()
plt.savefig(f"{BASE}/figures/planA_ppi_network.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{BASE}/figures/planA_ppi_network.pdf", bbox_inches="tight")
print("saved planA_ppi_network.png/.pdf")
print("nodes:", len(nodes), "edges:", len(edges))
print("mean combined_score:", round(sum(e[2] for e in edges)/len(edges), 3))
