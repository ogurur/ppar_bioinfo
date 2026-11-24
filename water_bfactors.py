# =====================================================================
# How firmly is each conserved water held?  (B-factors)
# ---------------------------------------------------------------------
# For every water site and every structure that has a water there (1 Å run),
# read the water's B-factor from the structure's mmCIF and normalise it:
#     B_norm = B(water) / mean B(Cα of the same LBD chain)
# B_norm < 1: the water is less mobile than the average protein backbone of
# its own structure; B_norm > 1: more mobile.  Normalising removes most of the
# difference between structures refined at different resolutions.
#
# Per site: median raw B, median B_norm, % of occurrences with raw B > 40 Å².
# Then: are highly conserved sites more ordered (lower B_norm) than rare ones?
#   Spearman ρ between site conservation (%) and site median B_norm,
#   Mann-Whitney U test high vs rare sites.
#
# Input (project folder):
#   conserved_align_waters_<subtype>.csv                (1 Å run)
#   data/results/water_sites_conservation_<subtype>.csv (plot_water_conservation.py)
#   data/raw/PPAR<subtype>/*.cif
# Output:
#   data/results/water_bfactors_<subtype>.csv        one row per water occurrence
#   data/results/water_site_bfactors_<subtype>.csv   one row per site
#   data/results/bfactor_stats.csv                   correlations per subtype
#   figures/bfactor_by_class.png
#   figures/bfactor_vs_conservation.png
#
# Run:  python water_bfactors.py        Needs: numpy, pandas, matplotlib, scipy
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR

HIGH, LOW = 70.0, 30.0       # conservation classes (% of structures)
B_FLAG = 40.0                # Å², poorly ordered water (project description)

SUBTYPES = {
    "alpha": {"label": "PPARα", "ref": "1k7l", "chain": "A", "raw": "data/raw/PPARalpha"},
    "gamma": {"label": "PPARγ", "ref": "1prg", "chain": "A", "raw": "data/raw/PPARgamma"},
    "delta": {"label": "PPARδ", "ref": "1gwx", "chain": "A", "raw": "data/raw/PPARdelta"},
}
WATER = {"HOH", "WAT", "DOD"}
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
      "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"}


def read_bfactors(path, chain):
    """{'ca_mean', 'ca_sd', 'waters': {resi: (B, occupancy)}} for one chain."""
    tags, rows, in_loop = [], [], False
    with open(path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                in_loop = True
                tags.append(s.split()[0][11:])
                continue
            if in_loop:
                if not s or s[0] in "#_" or s.startswith("loop_"):
                    if rows:
                        break
                    continue
                rows.append(s.split())
    c = {t: i for i, t in enumerate(tags)}
    resn_k = c.get("auth_comp_id", c["label_comp_id"])
    atom_k = c.get("auth_atom_id", c["label_atom_id"])
    ca, waters, m0 = [], {}, None
    for r in rows:
        if len(r) != len(tags) or r[c["auth_asym_id"]] != chain:
            continue
        if "pdbx_PDB_model_num" in c:
            m0 = m0 or r[c["pdbx_PDB_model_num"]]
            if r[c["pdbx_PDB_model_num"]] != m0:
                continue
        if "label_alt_id" in c and r[c["label_alt_id"]] not in (".", "?", "A", "1"):
            continue
        b = float(r[c["B_iso_or_equiv"]])
        if r[resn_k] in AA and r[atom_k] == "CA":
            ca.append(b)
        elif r[resn_k] in WATER and r[c["type_symbol"]] == "O":
            occ = float(r[c["occupancy"]]) if "occupancy" in c else 1.0
            waters.setdefault(r[c["auth_seq_id"]], (b, occ))
    ca = np.array(ca)
    return {"ca_mean": ca.mean() if len(ca) else np.nan,
            "ca_sd": ca.std() if len(ca) else np.nan, "waters": waters}


def klass(p):
    return "high" if p > HIGH else ("moderate" if p >= LOW else "rare")


def collect(key, cfg):
    sites = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                     f"water_sites_conservation_{key}.csv"),
                        dtype={"ref_water_resi": str})
    sites = sites[sites["n_1.0A"] > 0].copy()
    sites["class"] = sites["pct_1.0A"].map(klass)

    m = pd.read_csv(os.path.join(BASE_DIR, f"conserved_align_waters_{key}.csv"),
                    dtype={"water_resi": str, "chain": str})
    S = sites[["ref_x", "ref_y", "ref_z"]].to_numpy()
    d = np.linalg.norm(m[["x", "y", "z"]].to_numpy()[:, None] - S[None], axis=2)
    m["ref_water_resi"] = sites["ref_water_resi"].to_numpy()[d.argmin(1)]
    m["dist"] = d.min(1)
    m = m.sort_values("dist").drop_duplicates(["pdb_id", "ref_water_resi"])
    # the reference has a water at every one of its own sites
    ref = pd.DataFrame({"pdb_id": cfg["ref"], "chain": cfg["chain"],
                        "water_resi": sites["ref_water_resi"].values,
                        "ref_water_resi": sites["ref_water_resi"].values})
    m = pd.concat([ref, m[ref.columns]], ignore_index=True)

    rows, missing = [], 0
    for pdb_id, g in m.groupby("pdb_id"):
        path = os.path.join(BASE_DIR, cfg["raw"], f"{pdb_id}.cif")
        if not os.path.isfile(path):
            missing += len(g)
            continue
        info = {}
        for ch, gg in g.groupby("chain"):
            info = read_bfactors(path, ch)
            for _, r in gg.iterrows():
                if r["water_resi"] not in info["waters"]:
                    missing += 1
                    continue
                b, occ = info["waters"][r["water_resi"]]
                rows.append({"pdb_id": pdb_id, "chain": ch, "water_resi": r["water_resi"],
                             "ref_water_resi": r["ref_water_resi"], "B": b, "occupancy": occ,
                             "B_mean_CA": round(info["ca_mean"], 2),
                             "B_norm": round(b / info["ca_mean"], 3),
                             "B_z": round((b - info["ca_mean"]) / info["ca_sd"], 3)
                             if info["ca_sd"] > 0 else np.nan})
    occ = pd.DataFrame(rows)
    w_of = dict(zip(sites["ref_water_resi"], sites["W"]))
    occ.insert(0, "W", occ["ref_water_resi"].map(w_of))

    g = occ.groupby("ref_water_resi")
    per = pd.DataFrame({
        "n_occurrences": g.size(),
        "median_B": g["B"].median().round(1),
        "median_B_norm": g["B_norm"].median().round(3),
        "q1_B_norm": g["B_norm"].quantile(0.25).round(3),
        "q3_B_norm": g["B_norm"].quantile(0.75).round(3),
        "pct_B_over_40": (g["B"].apply(lambda b: (b > B_FLAG).mean()) * 100).round(1),
        "pct_occupancy_below_1": (g["occupancy"].apply(lambda o: (o < 1).mean()) * 100).round(1),
    }).reset_index()
    per = sites[["W", "ref_water_resi", "pct_1.0A", "class"]].merge(per, on="ref_water_resi")
    per = per.sort_values("pct_1.0A", ascending=False)
    return occ, per, missing


def test(per):
    rho, p = stats.spearmanr(per["pct_1.0A"], per["median_B_norm"])
    hi = per.loc[per["class"] == "high", "median_B_norm"]
    ra = per.loc[per["class"] == "rare", "median_B_norm"]
    u_p = stats.mannwhitneyu(hi, ra, alternative="two-sided").pvalue \
        if len(hi) >= 3 and len(ra) >= 3 else np.nan
    return rho, p, u_p, hi.median(), ra.median()


# ---------------------------------------------------------------------
INK, INK2, GRID, ACC, MUTED = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#a8a69f"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def plots(results, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    order = ["high", "moderate", "rare"]
    ymax = max(r[1]["median_B_norm"].max() for r in results.values()) * 1.08

    # 1 - site median B_norm by conservation class
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.3), dpi=300, sharey=True)
    for ax, (key, (occ, per, st)) in zip(axes, results.items()):
        data = [per.loc[per["class"] == c, "median_B_norm"].to_numpy() for c in order]
        ax.boxplot(data, positions=range(3), widths=0.5, showfliers=False,
                   medianprops=dict(color=INK, linewidth=1.5),
                   boxprops=dict(color=MUTED), whiskerprops=dict(color=MUTED),
                   capprops=dict(color=MUTED))
        for i, v in enumerate(data):
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14,
                       color=ACC if i == 0 else MUTED, alpha=0.85, zorder=3,
                       edgecolors="white", linewidths=0.4)
        ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.text(-0.45, 0.06, "dashed line: water as mobile as the average Cα",
                ha="left", va="bottom", fontsize=7, color=INK2)
        ax.set_xticks(range(3))
        ax.set_xticklabels([f"{c}\n(n={len(v)})" for c, v in zip(order, data)], fontsize=8)
        ax.set_ylim(0, ymax)
        rho, p, u_p, mh, mr = st
        ax.set_title(f"{SUBTYPES[key]['label']}", loc="left", fontsize=11, color=INK)
        ax.text(0.02, 0.98, f"high vs rare: {p_text(u_p)}" if not np.isnan(u_p)
                else "too few high sites to test", transform=ax.transAxes,
                fontsize=8, color=INK2, va="top")
        style(ax)
    axes[0].set_ylabel("Water B-factor / mean Cα B-factor\n(median per site)", fontsize=9, color=INK)
    fig.suptitle("Water mobility by conservation class (each dot = one water site, 1 Å)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(fig_dir, "bfactor_by_class.png"), facecolor="white")
    plt.close(fig)

    # 2 - conservation vs B_norm, one dot per site
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.3), dpi=300, sharey=True)
    for ax, (key, (occ, per, st)) in zip(axes, results.items()):
        col = np.where(per["class"] == "high", ACC, MUTED)
        ax.scatter(per["pct_1.0A"], per["median_B_norm"], s=16, c=col, alpha=0.85,
                   edgecolors="white", linewidths=0.4, zorder=3)
        hi = per[per["class"] == "high"]
        if len(hi) > 3:   # one label for a tight cluster
            ax.annotate(f"{hi['W'].iloc[0]}–{hi['W'].iloc[-1]}",
                        (hi["pct_1.0A"].mean(), hi["median_B_norm"].min()), xytext=(0, -14),
                        textcoords="offset points", ha="center", fontsize=7, color=INK)
        else:
            for i, (_, r) in enumerate(hi.sort_values("median_B_norm").iterrows()):
                ax.annotate(r["W"], (r["pct_1.0A"], r["median_B_norm"]),
                            xytext=(5, -7 if i == 0 else 5), textcoords="offset points",
                            va="center", fontsize=7, color=INK)
        ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.axvline(HIGH, color=GRID, linewidth=1)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, ymax)
        rho, p, *_ = st
        ax.set_title(f"{SUBTYPES[key]['label']}   Spearman ρ = {rho:.2f}, {p_text(p)}",
                     loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Site conservation (% of structures)", fontsize=9, color=INK)
        style(ax)
    axes[0].set_ylabel("Median water B / mean Cα B", fontsize=9, color=INK)
    fig.suptitle("More conserved water sites vs how mobile their waters are (blue = high, >70%)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(fig_dir, "bfactor_vs_conservation.png"), facecolor="white")
    plt.close(fig)


def main():
    res_dir = os.path.join(OUT_DIR, "data", "results")
    fig_dir = os.path.join(OUT_DIR, "figures")
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    results, srows = {}, []
    for key, cfg in SUBTYPES.items():
        occ, per, missing = collect(key, cfg)
        occ.to_csv(os.path.join(res_dir, f"water_bfactors_{key}.csv"), index=False)
        per.to_csv(os.path.join(res_dir, f"water_site_bfactors_{key}.csv"), index=False)
        st = test(per)
        results[key] = (occ, per, st)
        rho, p, u_p, mh, mr = st
        srows.append({"subtype": cfg["label"], "sites": len(per), "occurrences": len(occ),
                      "not_found": missing, "spearman_rho_pct_vs_Bnorm": round(rho, 3),
                      "spearman_p": p, "median_Bnorm_high": round(mh, 3),
                      "median_Bnorm_rare": round(mr, 3), "mannwhitney_p_high_vs_rare": u_p})
        print(f"== {cfg['label']}: {len(per)} sites, {len(occ)} waters, {missing} not found")
        print(f"   Spearman ρ(conservation, B_norm) = {rho:.2f} ({p_text(p)}); "
              f"median B_norm high {mh:.2f} vs rare {mr:.2f}"
              + (f" (Mann-Whitney {p_text(u_p)})" if not np.isnan(u_p) else ""))
        print(per[per["class"] != "rare"][["W", "pct_1.0A", "median_B", "median_B_norm",
                                           "pct_B_over_40"]].to_string(index=False))
    pd.DataFrame(srows).to_csv(os.path.join(res_dir, "bfactor_stats.csv"), index=False)
    plots(results, fig_dir)
    print("Tables in", res_dir, "| figures in", fig_dir)


if __name__ == "__main__":
    main()
