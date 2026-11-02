# =====================================================================
# Heatmap: conserved water sites x subtype (% of structures with a water)
# ---------------------------------------------------------------------
# Reads data/results/cross_subtype_occupancy.csv and
#       data/results/water_site_contacts_<subtype>.csv  (from water_contacts.py)
# Writes figures/cross_subtype_heatmap.png
# Run:  python plot_cross_subtype.py
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
MIN_ANY = 50.0      # show sites with >= this % in at least one subtype
SUB = [("alpha", "PPARα"), ("gamma", "PPARγ"), ("delta", "PPARδ")]

# one-hue sequential ramp (light = rarely present, dark = always present)
CMAP = LinearSegmentedColormap.from_list(
    "blue", ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
INK, INK_2 = "#0b0b0b", "#52514e"

res = os.path.join(BASE_DIR, "data", "results")
occ = pd.read_csv(os.path.join(res, "cross_subtype_occupancy.csv")).fillna("")
contacts = {}
for k, lab in SUB:
    c = pd.read_csv(os.path.join(res, f"water_site_contacts_{k}.csv")).fillna("")
    for _, r in c.iterrows():
        res_list = [x.rsplit(" ", 1)[0] for x in r["residues_contacted"].split("; ") if x]
        short = ", ".join(x[0] + x[1:3].lower() + x[3:] for x in res_list[:3]) or "waters only"
        contacts[f"{lab}-{r['W']}"] = (short, r["main_pattern"])

# merge duplicates: a site defined from two references gets both names
names = {}
for _, r in occ.iterrows():
    key = r["same_site_as"] or r["site"]
    names.setdefault(key, []).append(r["site"])
u = occ[occ["same_site_as"] == ""].copy()
cols = [f"pct_in_{lab}" for _, lab in SUB]
u = u[u[cols].max(axis=1) >= MIN_ANY]
u["n50"] = (u[cols] >= 50).sum(axis=1)
u = u.sort_values(["n50", "min_pct_all_subtypes"], ascending=False)

ORDER = {"PPARγ": 0, "PPARα": 1, "PPARδ": 2}
labels = []
for _, r in u.iterrows():
    group = sorted(names[r["site"]], key=lambda s: ORDER[s.split("-")[0]])
    alias = " = ".join(s.replace("PPAR", "") for s in group)
    short, pattern = contacts.get(group[0], ("", ""))   # residues of the first-named site
    labels.append(f"{alias}   ({short})")

M = u[cols].to_numpy(float)
fig_h = 0.32 * len(u) + 1.6
fig, ax = plt.subplots(figsize=(7.6, fig_h), dpi=300)
im = ax.imshow(M, cmap=CMAP, vmin=0, vmax=100, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                color="white" if v >= 55 else INK)
# separators between groups (3 / 2 / 1 subtypes >= 50 %)
n50 = u["n50"].to_numpy()
for i in range(1, len(n50)):
    if n50[i] != n50[i - 1]:
        ax.axhline(i - 0.5, color=INK, linewidth=1.2)
ax.set_xticks(range(len(SUB)))
ax.set_xticklabels([lab for _, lab in SUB], fontsize=9, color=INK)
ax.xaxis.tick_top()
ax.set_yticks(range(len(u)))
ax.set_yticklabels(labels, fontsize=7.5, color=INK)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label("% of structures with a water at the site (1 Å)", fontsize=8, color=INK_2)
cb.ax.tick_params(labelsize=7, colors=INK_2)
cb.outline.set_visible(False)
fig.suptitle("Conserved water sites across PPAR subtypes", x=0.02, ha="left",
             fontsize=11, color=INK)
ax.set_title("All structures superimposed on PPARγ 1PRG. Residues: numbering of the\n"
             "first-named site. Lines separate sites ≥50% in 3 / 2 / 1 subtypes.",
             fontsize=8, color=INK_2, loc="left", pad=24, x=-0.9)
fig.tight_layout()
out = os.path.join(OUT_DIR, "figures")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "cross_subtype_heatmap.png"), facecolor="white")
print("saved", os.path.join(out, "cross_subtype_heatmap.png"))
