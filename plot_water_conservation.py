# =====================================================================
# Water conservation plot: W1, W2, ... ranked from most to least conserved
# ---------------------------------------------------------------------
# Input  (from alignment_ppar_lbd*.py runs, in the project folder):
#   conserved_align_waters_<subtype>[_05|_2].csv   matched waters per cutoff
#   lbd_chains_PPAR<subtype>.csv                    -> total aligned structures
#   data/raw/PPAR<subtype>/<ref>.cif                -> reference water positions
#
# Each matched water is assigned to the NEAREST reference water.  A reference
# water = one site.  Conservation of a site = structures with a water there
# / all aligned structures (each structure counted once per site).
#
# Output:
#   figures/conservation_<subtype>_bars.png     W1..Wn, bars for each cutoff
#   figures/conservation_<subtype>_ranked.png   all sites, % vs rank
#   data/results/water_sites_conservation_<subtype>.csv   W label table
#   data/results/conservation_summary.csv       high/moderate/rare counts
#
# Run from the project folder:   python plot_water_conservation.py
# Needs: numpy, pandas, matplotlib
# =====================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------
BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR

MAIN_CUTOFF = 1.0        # W1, W2, ... are ranked by conservation at this cutoff
TOP_N_BARS = 40          # max waters shown in the bar chart (gamma has ~170)
HIGH, LOW = 0.70, 0.30   # manual's thresholds: high > 70 %, rare < 30 %

SUBTYPES = {
    "alpha": {
        "label": "PPARα", "ref_pdb": "1K7L", "ref_chain": "A",
        "ref_cif": "data/raw/PPARalpha/1k7l.cif",
        "chains_csv": "lbd_chains_PPARalpha.csv",
        "files": {0.5: "conserved_align_waters_alpha_05.csv",
                  1.0: "conserved_align_waters_alpha.csv",
                  2.0: "conserved_align_waters_alpha_2.csv"},
    },
    "gamma": {
        "label": "PPARγ", "ref_pdb": "1PRG", "ref_chain": "A",
        "ref_cif": "data/raw/PPARgamma/1prg.cif",
        "chains_csv": "lbd_chains_PPARgamma.csv",
        "files": {0.5: "conserved_align_waters_gamma_05.csv",
                  1.0: "conserved_align_waters_gamma.csv",
                  2.0: "conserved_align_waters_gamma_2.csv"},
    },
    "delta": {
        "label": "PPARδ", "ref_pdb": "1GWX", "ref_chain": "A",
        "ref_cif": "data/raw/PPARdelta/1gwx.cif",
        "chains_csv": "lbd_chains_PPARdelta.csv",
        "files": {0.5: "conserved_align_waters_delta_05.csv",
                  1.0: "conserved_align_waters_delta.csv",
                  2.0: "conserved_align_waters_delta_2.csv"},
    },
}

# one blue, light -> dark = small -> large cutoff (validated ordinal ramp)
CUTOFF_COLORS = {0.5: "#86b6ef", 1.0: "#2a78d6", 2.0: "#104281"}
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3de"


# ---------------------------------------------------------------------
# Reference waters straight from the mmCIF (same coordinates PyMOL used,
# since the reference never moves)
# ---------------------------------------------------------------------
def read_reference_waters(cif_path, chain):
    """{auth_seq_id: [(x, y, z), ...]} - several entries = alt. locations."""
    tags, rows, in_loop = [], [], False
    with open(cif_path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                in_loop = True
                tags.append(s.split()[0])
                continue
            if in_loop:
                if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                    if rows:
                        break
                    continue
                rows.append(s.split())
    col = {t: i for i, t in enumerate(tags)}

    def get(r, name):
        return r[col[f"_atom_site.{name}"]]

    waters = {}
    for r in rows:
        if len(r) != len(tags):
            continue
        if get(r, "label_comp_id") not in ("HOH", "WAT", "DOD"):
            continue
        if "pdbx_PDB_model_num" in " ".join(tags) and get(r, "pdbx_PDB_model_num") != "1":
            continue
        if get(r, "auth_asym_id") != chain:
            continue
        resi = get(r, "auth_seq_id")
        xyz = tuple(float(get(r, k)) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))
        waters.setdefault(resi, []).append(xyz)
    return waters


def n_aligned_structures(chains_csv):
    d = pd.read_csv(chains_csv, dtype=str).fillna("")
    return int((d["used_chains"].str.strip() != "").sum())


# ---------------------------------------------------------------------
# Site conservation for one subtype
# ---------------------------------------------------------------------
def site_table(cfg):
    ref = read_reference_waters(os.path.join(BASE_DIR, cfg["ref_cif"]), cfg["ref_chain"])
    resis = list(ref.keys())
    # every alt. location as a point, remembering which site it belongs to
    pts = np.array([xyz for r in resis for xyz in ref[r]])
    owner = np.array([r for r in resis for _ in ref[r]])
    total = n_aligned_structures(os.path.join(BASE_DIR, cfg["chains_csv"]))

    table = pd.DataFrame({
        "ref_water_resi": resis,
        "ref_x": [ref[r][0][0] for r in resis],
        "ref_y": [ref[r][0][1] for r in resis],
        "ref_z": [ref[r][0][2] for r in resis],
    })
    checks = []
    for cutoff, fname in sorted(cfg["files"].items()):
        w = pd.read_csv(os.path.join(BASE_DIR, fname))
        xyz = w[["x", "y", "z"]].to_numpy()
        d = np.linalg.norm(xyz[:, None, :] - pts[None, :, :], axis=2)
        w["ref_water_resi"] = owner[d.argmin(axis=1)]
        w["dist"] = d.min(axis=1)
        far = int((w["dist"] > cutoff + 0.02).sum())   # CSV coords are rounded to 0.01
        checks.append((cutoff, len(w), far))
        counts = w.groupby("ref_water_resi")["pdb_id"].nunique()
        n = table["ref_water_resi"].map(counts).fillna(0).astype(int)
        table[f"n_{cutoff}A"] = n
        table[f"pct_{cutoff}A"] = (100 * n / total).round(1)

    m = f"pct_{MAIN_CUTOFF}A"
    others = [c for c in table.columns if c.startswith("pct_") and c != m]
    table = table.sort_values([m] + sorted(others, reverse=True), ascending=False)
    table.insert(0, "W", [f"W{i}" for i in range(1, len(table) + 1)])
    table[f"class_{MAIN_CUTOFF}A"] = np.where(table[m] > 100 * HIGH, "high",
                                     np.where(table[m] >= 100 * LOW, "moderate", "rare"))
    table["n_total_structures"] = total
    return table.reset_index(drop=True), total, checks


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------
def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(0, 100)


def thresholds(ax, x_right):
    for y, txt in ((100 * HIGH, "high > 70%"), (100 * LOW, "rare < 30%")):
        ax.axhline(y, color=INK_2, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
        ax.text(x_right, y + 1.2, txt, ha="right", va="bottom", fontsize=8, color=INK_2)


def plot_bars(table, cfg, total, cutoffs, path):
    m = f"pct_{MAIN_CUTOFF}A"
    any_water = table[table[[f"n_{c}A" for c in cutoffs]].sum(axis=1) > 0]
    show = any_water.head(TOP_N_BARS)
    k, n = len(cutoffs), len(show)
    x = np.arange(n)
    width = 0.8 / k
    fig, ax = plt.subplots(figsize=(min(14, max(7, 0.3 * n + 2)), 4.6), dpi=300)
    for i, c in enumerate(cutoffs):
        ax.bar(x + (i - (k - 1) / 2) * width, show[f"pct_{c}A"], width * 0.9,
               color=CUTOFF_COLORS.get(c, "#2a78d6"), label=f"{c:g} Å",
               edgecolor="white", linewidth=0.4, zorder=2)
    style(ax)
    thresholds(ax, n - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(show["W"], rotation=90, fontsize=7.5, color=INK_2)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylabel("Structures with water at site (%)", color=INK, fontsize=9)
    ax.set_xlabel(f"Water site (ranked by conservation at {MAIN_CUTOFF:g} Å)",
                  color=INK, fontsize=9)
    extra = "" if n == len(any_water) else f", top {n} of {len(any_water)} sites shown"
    ax.set_title(f"{cfg['label']} water conservation, reference {cfg['ref_pdb']} "
                 f"(n = {total} structures{extra})",
                 loc="left", fontsize=11, color=INK)
    ax.legend(title="Water cutoff", fontsize=8, title_fontsize=8, frameon=False,
              ncol=k, loc="upper right", bbox_to_anchor=(1, 1.12))
    fig.tight_layout()
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def plot_ranked(table, cfg, total, cutoffs, path):
    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=300)
    for c in cutoffs:
        y = np.sort(table[f"pct_{c}A"].to_numpy())[::-1]
        y = y[y > 0]
        ax.step(np.arange(1, len(y) + 1), y, where="mid", linewidth=2,
                color=CUTOFF_COLORS.get(c, "#2a78d6"), label=f"{c:g} Å")
    style(ax)
    xmax = max(int((table[f"n_{c}A"] > 0).sum()) for c in cutoffs)
    thresholds(ax, xmax)
    ax.set_xlim(0.5, xmax + 0.5)
    ax.set_xlabel("Site rank (each cutoff ranked separately)", color=INK, fontsize=9)
    ax.set_ylabel("Structures with water at site (%)", color=INK, fontsize=9)
    ax.set_title(f"{cfg['label']}: all {len(table)} reference waters (n = {total} structures)",
                 loc="left", fontsize=11, color=INK)
    ax.legend(title="Water cutoff", fontsize=8, title_fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------
def main():
    fig_dir = os.path.join(OUT_DIR, "figures")
    res_dir = os.path.join(OUT_DIR, "data", "results")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)

    summary = []
    for name, cfg in SUBTYPES.items():
        table, total, checks = site_table(cfg)
        cutoffs = sorted(cfg["files"])
        print(f"\n== {cfg['label']}: {len(table)} reference waters, {total} aligned structures")
        for c, n_rows, far in checks:
            flag = "" if far == 0 else f"  [WARN] {far} waters farther than {c} A from any reference water"
            print(f"   {c:g} A: {n_rows} matched waters{flag}")

        table.to_csv(os.path.join(res_dir, f"water_sites_conservation_{name}.csv"), index=False)
        plot_bars(table, cfg, total, cutoffs, os.path.join(fig_dir, f"conservation_{name}_bars.png"))
        plot_ranked(table, cfg, total, cutoffs, os.path.join(fig_dir, f"conservation_{name}_ranked.png"))

        for c in cutoffs:
            p = table[f"pct_{c}A"]
            summary.append({
                "subtype": cfg["label"], "reference": cfg["ref_pdb"], "cutoff_A": c,
                "n_structures": total, "reference_waters": len(table),
                "sites_with_any_water": int((p > 0).sum()),
                "high_gt70": int((p > 100 * HIGH).sum()),
                "moderate_30_70": int(((p >= 100 * LOW) & (p <= 100 * HIGH)).sum()),
                "rare_lt30_nonzero": int(((p > 0) & (p < 100 * LOW)).sum()),
                "top_site": table.loc[p.idxmax(), "W"],
                "top_pct": float(p.max()),
            })
        top = table.head(5)[["W", "ref_water_resi"] + [f"pct_{c}A" for c in cutoffs]]
        print(top.to_string(index=False))

    pd.DataFrame(summary).to_csv(os.path.join(res_dir, "conservation_summary.csv"), index=False)
    print(f"\nFigures in {fig_dir}\nTables in {res_dir}")


if __name__ == "__main__":
    main()
