# -*- coding: utf-8 -*-
"""
Figure 8 — TF regulatory network of the beta-cell identity program.
PDX1/HNF1A -> identity/functional target genes.
Data: r3_inputs/tf_target_sets.csv (curated TF->target relations)
      ml_output/identity_log2fc_both.csv (real per-gene log2FC)
Node color = validation (GSE76894) log2FC (blue=down, red=up, gray=no probe).
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
import numpy as np

BASE = r"C:/Users/45902/WorkBuddy/2026-09-20-15-12-44"
OUT  = os.path.join(BASE, "figures")

# ---- load curated TF->target relations ----
sets = {}
for r in csv.DictReader(open(os.path.join(BASE, "r3_inputs", "tf_target_sets.csv"), encoding="utf-8")):
    sets.setdefault(r["set"], []).append(r["gene"])

pdx_targets  = sets["PDX1_targets"]        # PDX1 direct targets
hnf_targets  = sets["HNF1A_targets"]       # HNF1A direct targets
ident_net    = sets["beta_identity_network"]

# ---- load real log2FC (both cohorts) ----
fc = {}
for r in csv.DictReader(open(os.path.join(BASE, "ml_output", "identity_log2fc_both.csv"), encoding="utf-8")):
    fc[r["gene"]] = (float(r["log2FC_GSE76895"]), float(r["log2FC_GSE76894"]))

# validation log2FC used for node color (larger cohort, more consistent)
def val_fc(g):
    return fc[g][1] if g in fc else None

# ---- node universe ----
# Layer 1 TFs
tfs = ["HNF1A", "PDX1"]
# Layer 2 target genes: union of PDX1/HNF1A targets, keep only meaningful
#   (exclude pure-TF nodes that are themselves layer-1 to avoid clutter; keep HNF4A as a key target)
targets = sorted(set(pdx_targets) | set(hnf_targets))
# PDX1 & HNF1A themselves appear in target lists (HNF1A->PDX1 edge; keep PDX1 as TF node)
targets = [g for g in targets if g not in tfs]

# distinguish which TF regulates each target
pdx_set = set(pdx_targets)
hnf_set = set(hnf_targets)

# ---- manual layout (x, y) ----
# TFs on left; targets on right, grouped
pos = {
    "HNF1A": (0.0, 0.55),
    "PDX1":  (0.0, -0.55),
}
# key targets with real log2FC placed right
target_pos = {
    "SLC2A2":  (1.05, 0.80),   # hub (both TFs), strongest down
    "GCK":     (1.05, 0.40),
    "INS":     (1.05, 0.00),
    "IAPP":    (1.05, -0.40),
    "NKX6-1":  (1.05, -0.80),
    "HNF4A":   (2.10, 0.70),
    "SLC30A8": (2.10, 0.30),
    "ABCC8":   (2.10, -0.10),
    "KCNJ11":  (2.10, -0.50),
    "PKLR":    (2.10, -0.90),
}
# column 3: remaining targets, evenly spaced (HNF4G, FOXA3, HNF1B, MAFA, PAX4, PC, G6PC2, GLUT2)
col3 = ["HNF4G", "FOXA3", "HNF1B", "MAFA", "PAX4", "PC", "G6PC2", "GLUT2"]
for i, g in enumerate(col3):
    target_pos[g] = (3.15, 0.95 - 0.28*i)
for g in targets:
    if g not in target_pos:
        target_pos[g] = (3.15, 0.0)  # fallback (shouldn't happen)
# place remaining (HNF4G, FOXA3, HNF1B) in a third sub-column, top
extra = [g for g in targets if g not in target_pos]
for i, g in enumerate(extra):
    target_pos[g] = (3.15, 0.95 - 0.42*i)

# ---- colors ----
# validation log2FC colormap: down = blue, up = red, no-probe = gray
from matplotlib.colors import LinearSegmentedColormap
def node_color(g):
    v = val_fc(g)
    if v is None:
        return "#bdbdbd"   # gray = no GPL570 probe
    # clamp to [-0.5, 0.1]; down negative
    t = max(-0.5, min(0.1, v))
    # RdBu colormap: 0.0 = deep red, 0.5 = white, 1.0 = deep blue
    # t = -0.5 (strongest down) -> 1.0 (deep blue); t = 0 -> 0.5 (white); t = +0.1 -> 0.4 (light red)
    u = 0.5 - t
    return plt.cm.RdBu(min(1.0, max(0.0, u)))

fig, ax = plt.subplots(figsize=(11.5, 6.2))

# ---- edges ----
edge_pdx_color = "#2c7fb8"   # PDX1 -> target (blue)
edge_hnf_color = "#e08214"   # HNF1A -> target (orange)
edge_pdx = set(pdx_targets)
edge_hnf = set(hnf_targets)

# HNF1A -> PDX1 (TF-TF cascade, thick)
ax.add_patch(FancyArrowPatch(pos["HNF1A"], pos["PDX1"],
    arrowstyle="-|>", mutation_scale=18, lw=2.6,
    color="#5a5a5a", alpha=0.9, connectionstyle="arc3,rad=-0.25",
    linestyle="--", zorder=1))

# TF -> target edges
for g in targets:
    gp = target_pos[g]
    if g in edge_pdx:
        ax.add_patch(FancyArrowPatch(pos["PDX1"], gp,
            arrowstyle="-|>", mutation_scale=13, lw=1.5,
            color=edge_pdx_color, alpha=0.55, zorder=1,
            connectionstyle="arc3,rad=0.0"))
    if g in edge_hnf:
        ax.add_patch(FancyArrowPatch(pos["HNF1A"], gp,
            arrowstyle="-|>", mutation_scale=13, lw=1.5,
            color=edge_hnf_color, alpha=0.55, zorder=1,
            connectionstyle="arc3,rad=0.0"))

# ---- nodes ----
def draw_node(g, p, r=0.13, is_tf=False):
    x, y = p
    c = node_color(g)
    if is_tf:
        # TFs: rounded rectangle
        ax.add_patch(FancyBboxPatch((x-0.30, y-0.13), 0.60, 0.26,
            boxstyle="round,pad=0.02", fc=c, ec="black", lw=1.5, zorder=3))
        # adaptive text color for readability
        rgb = matplotlib.colors.to_rgb(c)
        lum = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
        tcol = "white" if lum < 0.5 else "black"
        ax.text(x, y, g, ha="center", va="center", fontsize=11,
                fontweight="bold", color=tcol, zorder=4)
    else:
        ax.add_patch(Circle((x, y), r, fc=c, ec="black", lw=1.2, zorder=3))
        ax.text(x, y, g, ha="center", va="center", fontsize=8,
                color="black", zorder=4, fontweight="bold")

draw_node("HNF1A", pos["HNF1A"], is_tf=True)
draw_node("PDX1",  pos["PDX1"],  is_tf=True)
for g in targets:
    draw_node(g, target_pos[g])

# gold highlight ring on SLC2A2 (the only shared target of both TFs)
sx, sy = target_pos["SLC2A2"]
ax.add_patch(Circle((sx, sy), 0.175, fc="none", ec="#d4a017", lw=2.2, zorder=2.5))
ax.annotate("shared hub:\ntargeted by both TFs", xy=(sx+0.17, sy+0.10),
            xytext=(sx+0.52, sy+0.28), fontsize=7.5, style="italic",
            color="#8a6d1a", ha="left",
            arrowprops=dict(arrowstyle="-", color="#d4a017", lw=1.0))

# ---- legend / annotations ----
# colorbar for log2FC (RdBu_r + TwoSlopeNorm centered at 0:
#   -0.5 = deep blue (strongest down), 0 = white, +0.1 = deep red (up))
from matplotlib.colors import TwoSlopeNorm
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02, shrink=0.75)
cbar.set_label("log$_2$FC (T2DM vs ND, GSE76894)", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# legend for edges
from matplotlib.lines import Line2D
handles = [
    Line2D([0],[0], color=edge_pdx_color, lw=2, label="PDX1 target"),
    Line2D([0],[0], color=edge_hnf_color, lw=2, label="HNF1A target"),
    Line2D([0],[0], color="#5a5a5a", lw=2, ls="--", label="TF\u2192TF cascade"),
    Line2D([0],[0], marker="o", color="none", markerfacecolor="#bdbdbd",
           markeredgecolor="black", markersize=10, label="not on common probe panel"),
]
ax.legend(handles=handles, loc="lower left", fontsize=8.5, frameon=True, framealpha=0.9)

ax.text(2.3, 1.42, "Blue = downregulated in T2DM (bulk, GSE76894); gold ring = targeted by both TFs",
        fontsize=8.5, color="#333333", style="italic", ha="center")
ax.set_title("Transcriptional cascade of the \u03b2-cell identity program\n"
             "PDX1 / HNF1A target genes collectively downregulated in T2DM",
             fontsize=12, fontweight="bold", pad=8)

ax.set_xlim(-0.5, 4.1)
ax.set_ylim(-1.45, 1.52)
ax.axis("off")

out_pdf = os.path.join(OUT, "planA_tf_network.pdf")
out_png = os.path.join(OUT, "planA_tf_network.png")
plt.tight_layout()
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print("saved:", out_pdf)
print("saved:", out_png)

# ---- also dump a node list for Supplementary Table ----
out_csv = os.path.join(BASE, "ml_output", "tf_network_nodes.csv")
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["node", "type", "regulator(s)", "log2FC_GSE76895", "log2FC_GSE76894", "in_identity_network"])
    for g in ["HNF1A", "PDX1"]:
        w.writerow([g, "TF", "", (fc[g][0] if g in fc else ""), (fc[g][1] if g in fc else ""), g in ident_net])
    for g in targets:
        regs = []
        if g in edge_pdx: regs.append("PDX1")
        if g in edge_hnf: regs.append("HNF1A")
        v = fc[g] if g in fc else ("", "")
        w.writerow([g, "target", "|".join(regs), v[0], v[1], g in ident_net])
print("saved:", out_csv)
