# =====================================================================
# Structures with no water at any reference site: why?
# ---------------------------------------------------------------------
# "Zero-match" = an aligned structure (reference excluded) with no water
# within 1 Å of ANY reference water of its subtype.
# For these structures: resolution, number of waters modelled in the LBD,
# space group and depositing group, compared with all other structures.
# Test: are zero-match structures at worse resolution? (Mann-Whitney U)
#
# Input (project folder):
#   data/results/resolution_per_structure_<subtype>.csv (resolution_vs_conservation.py)
#   data/metadata.csv
# Output:
#   data/results/zero_match_structures.csv     the zero-match list with details
#   data/results/zero_match_stats.csv          counts and tests per subtype
#   figures/zero_match_resolution.png
#
# Run:  python zero_match_structures.py      Needs: numpy, pandas, matplotlib, scipy
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR
SUB = [("alpha", "PPARα"), ("gamma", "PPARγ"), ("delta", "PPARδ")]
INK, INK2, GRID, ACC, MUTED = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#a8a69f"


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def main():
    res = os.path.join(BASE_DIR, "data", "results")
    meta = pd.read_csv(os.path.join(BASE_DIR, "data", "metadata.csv"))
    meta["pdb_id"] = meta["pdb_id"].str.lower()
    meta = meta.set_index("pdb_id")

    data, zero_rows, stat_rows = {}, [], []
    for key, lab in SUB:
        d = pd.read_csv(os.path.join(res, f"resolution_per_structure_{key}.csv"))
        d["zero"] = ~d["any_site_occupied"].astype(bool)
        data[key] = d
        z, o = d[d["zero"]], d[~d["zero"]]
        u_p = stats.mannwhitneyu(z["resolution_A"], o["resolution_A"],
                                 alternative="greater").pvalue if len(z) >= 3 else np.nan
        no_w = z["n_waters_in_lbd"].fillna(-1)
        stat_rows.append({
            "subtype": lab, "structures": len(d), "zero_match": len(z),
            "zero_match_pct": round(100 * len(z) / len(d), 1),
            "median_res_zero": z["resolution_A"].median(),
            "median_res_others": o["resolution_A"].median(),
            "zero_match_worse_than_2.5A": int((z["resolution_A"] > 2.5).sum()),
            "zero_match_with_no_modelled_waters": int((no_w == 0).sum()),
            "zero_match_waters_unknown": int((no_w < 0).sum()),
            "zero_match_with_some_waters": int((no_w > 0).sum()),
            "mannwhitney_p_zero_worse": u_p,
        })
        for _, r in z.iterrows():
            m = meta.loc[r.pdb_id] if r.pdb_id in meta.index else {}
            zero_rows.append({
                "subtype": lab, "pdb_id": r.pdb_id, "resolution_A": r.resolution_A,
                "n_waters_in_lbd": r.n_waters_in_lbd,
                "space_group": m.get("space_group", ""), "depositing_group": m.get("depositing_group", ""),
                "deposit_date": m.get("deposit_date", ""), "ligands_present": m.get("ligands_present", ""),
            })
        print(f"== {lab}: {len(z)} of {len(d)} structures have no water at any site; "
              f"median resolution {z['resolution_A'].median():.2f} vs {o['resolution_A'].median():.2f} Å "
              f"({p_text(u_p)})")

    zdf = pd.DataFrame(zero_rows).sort_values(["subtype", "resolution_A"])
    zdf.to_csv(os.path.join(res, "zero_match_structures.csv"), index=False)
    st = pd.DataFrame(stat_rows)
    st.to_csv(os.path.join(res, "zero_match_stats.csv"), index=False)
    print(st.drop(columns="mannwhitney_p_zero_worse").to_string(index=False))
    print(zdf[["subtype", "pdb_id", "resolution_A", "n_waters_in_lbd", "depositing_group"]]
          .to_string(index=False))

    # ---- figure: resolution of zero-match vs other structures -------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(3, 1, figsize=(9, 6.2), dpi=300, sharex=True)
    for ax, (key, lab), srow in zip(axes, SUB, stat_rows):
        d = data[key]
        for y, grp, col in ((1, d[~d["zero"]], MUTED), (0, d[d["zero"]], ACC)):
            yy = y + rng.uniform(-0.18, 0.18, len(grp))
            nw = grp["n_waters_in_lbd"].fillna(-1).to_numpy()
            filled = nw != 0
            ax.scatter(grp["resolution_A"][filled], yy[filled], s=16, color=col, alpha=0.75,
                       edgecolors="white", linewidths=0.4, zorder=3)
            ax.scatter(grp["resolution_A"][~filled], yy[~filled], s=18, facecolors="none",
                       edgecolors=col, linewidths=1.1, zorder=3)
            ax.plot([grp["resolution_A"].median()] * 2, [y - 0.3, y + 0.3], color=INK,
                    linewidth=2, zorder=4)
        ax.axvline(2.5, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.set_yticks([0, 1])
        ax.set_yticklabels([f"no water at any site (n={srow['zero_match']})",
                            f"water at ≥1 site (n={srow['structures'] - srow['zero_match']})"],
                           fontsize=8, color=INK)
        ax.set_ylim(-0.55, 1.55)
        ax.set_title(f"{lab}   zero-match structures at worse resolution: "
                     f"{p_text(srow['mannwhitney_p_zero_worse'])}",
                     loc="left", fontsize=9.5, color=INK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(INK2)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
    axes[-1].set_xlabel("Resolution (Å, lower = better)", fontsize=9, color=INK)
    axes[0].text(2.52, 1.5, "2.5 Å", fontsize=7, color=INK2, va="top")
    fig.suptitle("Structures with no water at any reference site are mostly low-resolution\n"
                 "dot = structure · open circle = no waters modelled in the LBD at all · "
                 "black bar = median", x=0.01, ha="left", fontsize=10.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(OUT_DIR, "figures")
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "zero_match_resolution.png"), facecolor="white")
    print("Tables in", res, "| figure in", out)


if __name__ == "__main__":
    main()
