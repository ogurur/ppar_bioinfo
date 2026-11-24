# =====================================================================
# Does water conservation depend on crystallographic resolution?
# ---------------------------------------------------------------------
# For every aligned structure (reference excluded):
#   - its resolution (data/metadata.csv);
#   - the % of the subtype's conserved sites (>= MIN_PCT at 1 Å) where it
#     has a water  -> "conserved sites occupied";
#   - for comparison, the total number of waters modelled in its LBD.
# Spearman ρ (with p) between resolution and each measure, per subtype.
# If conserved sites are real, they should be found at most resolutions:
# a weak or no correlation, unlike the total water count, which falls as
# resolution gets worse.
#
# Also, for each conserved site: % of structures with a water there in three
# resolution groups (≤2.0, 2.0–2.5, >2.5 Å).
#
# Input (project folder):
#   data/metadata.csv, lbd_chains_PPAR<subtype>.csv,
#   conserved_align_waters_<subtype>.csv (1 Å run),
#   data/results/water_sites_conservation_<subtype>.csv
# Output:
#   data/results/resolution_per_structure_<subtype>.csv
#   data/results/resolution_per_site_<subtype>.csv
#   data/results/resolution_stats.csv
#   figures/resolution_vs_conservation.png
#
# Run:  python resolution_vs_conservation.py   Needs: numpy, pandas, matplotlib, scipy
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR

MIN_PCT = 30.0                     # conserved site = present in >= 30% at 1 Å
BINS = [0, 2.0, 2.5, 10]           # resolution groups for the per-site table
BIN_LABELS = ["≤2.0 Å", "2.0–2.5 Å", ">2.5 Å"]

SUBTYPES = {
    "alpha": {"label": "PPARα", "ref": "1k7l", "meta": "PPARalpha"},
    "gamma": {"label": "PPARγ", "ref": "1prg", "meta": "PPARgamma"},
    "delta": {"label": "PPARδ", "ref": "1gwx", "meta": "PPARdelta"},
}


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def per_structure(key, cfg, meta):
    sites = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                     f"water_sites_conservation_{key}.csv"),
                        dtype={"ref_water_resi": str})
    occupied_sites = sites[sites["n_1.0A"] > 0]
    cons = sites[sites["pct_1.0A"] >= MIN_PCT]

    m = pd.read_csv(os.path.join(BASE_DIR, f"conserved_align_waters_{key}.csv"),
                    dtype={"water_resi": str})
    S = occupied_sites[["ref_x", "ref_y", "ref_z"]].to_numpy()
    d = np.linalg.norm(m[["x", "y", "z"]].to_numpy()[:, None] - S[None], axis=2)
    m["ref_water_resi"] = occupied_sites["ref_water_resi"].to_numpy()[d.argmin(1)]
    hits = m.drop_duplicates(["pdb_id", "ref_water_resi"])

    chains = pd.read_csv(os.path.join(BASE_DIR, f"lbd_chains_PPAR{key}.csv"), dtype=str).fillna("")
    ids = chains.loc[chains["used_chains"] != "", "pdb_id"].str.lower()
    ids = ids[ids != cfg["ref"]]

    rows = []
    for pid in ids:
        h = hits[hits["pdb_id"] == pid]
        n_cons = h["ref_water_resi"].isin(cons["ref_water_resi"]).sum()
        md = meta.loc[pid.upper()]
        rows.append({"pdb_id": pid, "resolution_A": md["resolution_A"],
                     "n_waters_in_lbd": md["n_waters_in_lbd"],
                     "conserved_sites_occupied": int(n_cons),
                     "pct_conserved_sites_occupied": round(100 * n_cons / len(cons), 1),
                     "any_site_occupied": len(h) > 0})
    return pd.DataFrame(rows), cons, hits


def per_site(cons, hits, ps):
    ps = ps.copy()
    ps["res_group"] = pd.cut(ps["resolution_A"], BINS, labels=BIN_LABELS)
    groups = ps.groupby("res_group", observed=False)["pdb_id"].apply(set)
    rows = []
    for _, s in cons.iterrows():
        have = set(hits.loc[hits["ref_water_resi"] == s.ref_water_resi, "pdb_id"])
        row = {"W": s.W, "ref_water_resi": s.ref_water_resi, "pct_all": s["pct_1.0A"]}
        for lab in BIN_LABELS:
            g = groups.get(lab, set())
            row[f"pct_{lab}"] = round(100 * len(g & have) / len(g), 1) if g else np.nan
            row[f"n_{lab}"] = len(g)
        # within-site test: is presence related to resolution?
        x = ps["resolution_A"].to_numpy()
        y = ps["pdb_id"].isin(have).to_numpy().astype(int)
        rho, p = stats.spearmanr(x, y) if y.std() > 0 else (np.nan, np.nan)
        row["spearman_rho_presence_vs_res"] = round(rho, 3)
        row["spearman_p"] = p
        rows.append(row)
    return pd.DataFrame(rows)


INK, INK2, GRID, ACC, MUTED = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#a8a69f"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def trend(ax, x, y, color, step=0.25):
    """Median (dot) and middle 50% (bar) of y in resolution groups of `step` Å."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    edges = np.arange(np.floor(x.min() / step) * step, x.max() + step, step)
    for lo, hi in zip(edges[:-1], edges[1:]):
        v = y[(x >= lo) & (x < hi)]
        if len(v) < 3:
            continue
        mid = (lo + hi) / 2
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([mid, mid], [q1, q3], color=color, linewidth=2, zorder=4)
        ax.scatter([mid], [med], s=30, color=color, zorder=5, edgecolors="white", linewidths=0.8)


def jitter(x, rng, w=0.02):
    """Small horizontal spread so structures with the same resolution do not stack."""
    return np.asarray(x, float) + rng.uniform(-w, w, len(x))


def plot(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.4), dpi=300, sharex=True)
    rng = np.random.default_rng(0)
    for j, (key, (ps, st)) in enumerate(results.items()):
        lab = SUBTYPES[key]["label"]
        # top: conserved sites occupied
        ax = axes[0, j]
        ax.scatter(jitter(ps["resolution_A"], rng), ps["pct_conserved_sites_occupied"], s=14, color=ACC,
                   alpha=0.55, edgecolors="white", linewidths=0.3, zorder=3)
        trend(ax, ps["resolution_A"], ps["pct_conserved_sites_occupied"], INK)
        ax.set_ylim(-3, 103)
        ax.set_title(f"{lab}   ρ = {st['rho_cons']:.2f}, {p_text(st['p_cons'])}",
                     loc="left", fontsize=10, color=INK)
        style(ax)
        # bottom: all modelled waters, for comparison
        ax = axes[1, j]
        ok = ps["n_waters_in_lbd"].notna()
        ax.scatter(jitter(ps.loc[ok, "resolution_A"], rng), ps.loc[ok, "n_waters_in_lbd"], s=14, color=MUTED,
                   alpha=0.6, edgecolors="white", linewidths=0.3, zorder=3)
        trend(ax, ps.loc[ok, "resolution_A"], ps.loc[ok, "n_waters_in_lbd"], INK)
        ax.set_title(f"{lab}   ρ = {st['rho_all']:.2f}, {p_text(st['p_all'])}",
                     loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Resolution (Å, lower = better)", fontsize=9, color=INK)
        style(ax)
    axes[0, 0].set_ylabel(f"Conserved sites (≥{MIN_PCT:g}%) with a water\n(% of the subtype's sites)",
                          fontsize=9, color=INK)
    axes[1, 0].set_ylabel("All waters modelled in the LBD\n(count)", fontsize=9, color=INK)
    fig.suptitle("Resolution vs conserved waters (top) and vs all modelled waters (bottom); "
                 "one dot per structure;\nblack = median and middle 50% of each 0.25 Å resolution group",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def main():
    res_dir = os.path.join(OUT_DIR, "data", "results")
    fig_dir = os.path.join(OUT_DIR, "figures")
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    meta = pd.read_csv(os.path.join(BASE_DIR, "data", "metadata.csv")).set_index("pdb_id")
    results, srows = {}, []
    for key, cfg in SUBTYPES.items():
        ps, cons, hits = per_structure(key, cfg, meta)
        site = per_site(cons, hits, ps)
        ps.to_csv(os.path.join(res_dir, f"resolution_per_structure_{key}.csv"), index=False)
        site.to_csv(os.path.join(res_dir, f"resolution_per_site_{key}.csv"), index=False)
        rc, pc = stats.spearmanr(ps["resolution_A"], ps["pct_conserved_sites_occupied"])
        ok = ps["n_waters_in_lbd"].notna()
        ra, pa = stats.spearmanr(ps.loc[ok, "resolution_A"], ps.loc[ok, "n_waters_in_lbd"])
        st = {"rho_cons": rc, "p_cons": pc, "rho_all": ra, "p_all": pa}
        results[key] = (ps, st)
        n_sig = int((site["spearman_p"] < 0.05).sum())
        srows.append({"subtype": cfg["label"], "structures": len(ps), "conserved_sites": len(cons),
                      "spearman_rho_res_vs_conserved_occupied": round(rc, 3), "p": pc,
                      "spearman_rho_res_vs_all_waters": round(ra, 3), "p_all_waters": pa,
                      "sites_with_presence_related_to_resolution_p<0.05": n_sig})
        print(f"== {cfg['label']}: {len(ps)} structures, {len(cons)} conserved sites")
        print(f"   resolution vs conserved sites occupied: ρ = {rc:.2f} ({p_text(pc)})")
        print(f"   resolution vs all modelled waters:      ρ = {ra:.2f} ({p_text(pa)})")
        print(f"   sites whose presence depends on resolution (p<0.05): {n_sig} of {len(site)}")
        print(site.head(8)[["W", "pct_all"] + [f"pct_{b}" for b in BIN_LABELS]
                           + ["spearman_rho_presence_vs_res", "spearman_p"]]
              .round(3).to_string(index=False))
    pd.DataFrame(srows).to_csv(os.path.join(res_dir, "resolution_stats.csv"), index=False)
    plot(results, os.path.join(fig_dir, "resolution_vs_conservation.png"))
    print("Tables in", res_dir, "| figure in", fig_dir)


if __name__ == "__main__":
    main()
