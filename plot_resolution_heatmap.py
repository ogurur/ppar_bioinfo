# =====================================================================
# Heatmap: conserved water sites x resolution group
# ---------------------------------------------------------------------
# Rows = conserved sites (>= 30% at 1 Å), ranked; columns = all structures,
# then resolution groups ≤2.0, 2.0–2.5, >2.5 Å.  Cell = % of structures in
# that group with a water at the site.  Row labels carry the nearest helix.
#
# Input:  data/results/resolution_per_site_<subtype>.csv  (resolution_vs_conservation.py)
#         data/results/water_site_helices_<subtype>.csv   (label_helices.py, optional)
# Output: figures/resolution_heatmap.png
# Run:    python plot_resolution_heatmap.py
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR
SUB = [("alpha", "PPARα"), ("gamma", "PPARγ"), ("delta", "PPARδ")]
GROUPS = ["≤2.0 Å", "2.0–2.5 Å", ">2.5 Å"]
CMAP = LinearSegmentedColormap.from_list(
    "blue", ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
INK, INK2 = "#0b0b0b", "#52514e"

res = os.path.join(BASE_DIR, "data", "results")
tables = {}
for key, lab in SUB:
    t = pd.read_csv(os.path.join(res, f"resolution_per_site_{key}.csv"),
                    dtype={"ref_water_resi": str})
    hx = os.path.join(res, f"water_site_helices_{key}.csv")
    if os.path.isfile(hx):
        h = pd.read_csv(hx)[["W", "nearest_helix"]]
        t = t.merge(h, on="W", how="left")
    else:
        t["nearest_helix"] = ""
    tables[key] = t.sort_values("pct_all", ascending=False).reset_index(drop=True)

nmax = max(len(t) for t in tables.values())
fig = plt.figure(figsize=(11, 0.28 * nmax + 2.2), dpi=300)
gs = fig.add_gridspec(nmax, 3, wspace=0.55, left=0.07, right=0.9, top=1 - 1.25 / (0.28 * nmax + 2.2),
                      bottom=0.04)
for j, (key, lab) in enumerate(SUB):
    t = tables[key]
    ax = fig.add_subplot(gs[:len(t), j])
    cols = ["pct_all"] + [f"pct_{g}" for g in GROUPS]
    M = t[cols].to_numpy(float)
    im = ax.imshow(M, cmap=CMAP, vmin=0, vmax=100, aspect="auto")
    for i in range(M.shape[0]):
        for k in range(M.shape[1]):
            v = M[i, k]
            if np.isnan(v):
                continue
            ax.text(k, i, f"{v:.0f}", ha="center", va="center", fontsize=7,
                    color="white" if v >= 55 else INK)
    ax.axvline(0.5, color="white", linewidth=4)          # separate "all" from the groups
    n = [t[f"n_{g}"].iloc[0] for g in GROUPS]
    ax.set_xticks(range(4))
    ax.set_xticklabels(["all"] + [f"{g}\nn={k}" for g, k in zip(GROUPS, n)], fontsize=7.5, color=INK)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{w}  {h}" if isinstance(h, str) and h else w
                        for w, h in zip(t["W"], t["nearest_helix"])], fontsize=7.5, color=INK)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(lab, loc="left", fontsize=11, color=INK, pad=26)
cax = fig.add_axes([0.925, 0.55, 0.012, 0.3])
cb = fig.colorbar(im, cax=cax)
cb.set_label("% of structures with a water at the site (1 Å)", fontsize=8, color=INK2)
cb.ax.tick_params(labelsize=7, colors=INK2)
cb.outline.set_visible(False)
fig.suptitle("Conserved water sites by crystallographic resolution "
             "(rows = sites ≥30%, ranked; label = nearest helix)",
             x=0.01, ha="left", fontsize=12, color=INK)
out = os.path.join(OUT_DIR, "figures")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "resolution_heatmap.png"), facecolor="white")
print("saved", os.path.join(out, "resolution_heatmap.png"))
