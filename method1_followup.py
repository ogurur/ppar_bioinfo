# =====================================================================
# Supervisor's tasks 1-3, redone on the Method 1 clusters
# ---------------------------------------------------------------------
# Same questions as label_helices.py, water_bfactors.py and
# resolution_vs_conservation.py, but a "site" is now a Method 1 cluster
# (method1_cluster.py, greedy 1.5 Å, core-helix alignment) instead of a
# reference water, so every analysis uses one definition of a site.
#
#   Task 1  which helix / loop is each conserved cluster next to?
#           (residues within CONTACT Å of the cluster centre in the reference)
#   Task 2  are conserved clusters more ordered?  Median B_norm per cluster
#           (water B / mean Cα B of its chain) vs conservation, Spearman +
#           Mann-Whitney high vs rare.
#   Task 3  does conservation depend on resolution?  Per structure: % of the
#           conserved clusters (>= 30 %) where it has a water, vs resolution;
#           per cluster: occupancy in ≤2.0 / 2.0-2.5 / >2.5 Å structures.
#
# Input:  data/processed/method1/*, data/results/method1/clusters_<subtype>.csv,
#         data/raw/..., method1_cluster.py and label_helices.py (imported)
# Output: data/results/method1/followup_helices_<subtype>.csv
#         data/results/method1/followup_bfactor_<subtype>.csv
#         data/results/method1/followup_resolution_<subtype>.csv   (per structure)
#         data/results/method1/followup_resolution_sites_<subtype>.csv
#         data/results/method1/followup_stats.csv
#         figures/m1_sites_per_helix.png, m1_bfactor_vs_conservation.png,
#         figures/m1_resolution_vs_conservation.png, m1_resolution_heatmap.png
#
# Run:  python method1_followup.py      (after method1_prepare.py and method1_cluster.py)
# =====================================================================
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
sys.argv = sys.argv[:1] + [BASE_DIR]          # the imported scripts read argv[1]
import method1_cluster as m1                   # noqa: E402
import label_helices as lh                     # noqa: E402

RES = os.path.join(BASE_DIR, "data", "results", "method1")
FIG = os.path.join(BASE_DIR, "figures")
CONTACT = 4.0
MIN_STRUCT = 3            # task 2: clusters seen in at least this many structures
RES_BINS = [0, 2.0, 2.5, 10]
RES_LABELS = ["≤2.0 Å", "2.0–2.5 Å", ">2.5 Å"]
HEAT_N = 25               # clusters shown per subtype in the heatmap
HIGH, LOW = m1.HIGH, m1.LOW
SUB = {"alpha": ("PPARα", "1k7l", "data/raw/PPARalpha"),
       "gamma": ("PPARγ", "1prg", "data/raw/PPARgamma"),
       "delta": ("PPARδ", "1gwx", "data/raw/PPARdelta")}
INK, INK2, GRID, ACC, MUTED = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#a8a69f"


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def membership(key, ref):
    """Re-run the main clustering and return waters with their cluster label,
    the full cluster table, and the atlas names (C1, C2, ...)."""
    W, T = m1.load(key)
    X, ok, rmsd = m1.place(W, T, ref, "core")
    good = W.pdb_id.map(lambda p: rmsd.get(p, np.inf) <= m1.MAX_RMSD).to_numpy() & ok
    Wm, Xm = W[good].reset_index(drop=True), X[good]
    labels = m1.cluster(Xm, *m1.MAIN)
    Wm["label"] = labels
    n_total = Wm.pdb_id.nunique()
    cl = m1.summarise(labels, Xm, Wm.pdb_id.to_numpy(), n_total)
    atlas = pd.read_csv(os.path.join(RES, f"clusters_{key}.csv"))
    # name clusters as in the atlas by matching centres
    A = atlas[["cx", "cy", "cz"]].to_numpy()
    C = cl[["cx", "cy", "cz"]].to_numpy()
    d = np.linalg.norm(C[:, None] - A[None], axis=2)
    cl["C"] = np.where(d.min(1) < 0.05, atlas["C"].to_numpy()[d.argmin(1)], "")
    cl = cl.merge(atlas[["C", "pct_nr", "nearest_ref_site", "in_reference"]], on="C", how="left")
    return Wm, cl, n_total


# ---------------------------------------------------------------- task 1
def helices(key, raw, ref, cons, res_g):
    atoms, res = lh.read_chain(os.path.join(BASE_DIR, raw, f"{ref}.cif"), "A")
    gmap = {r[0]: r[0] for r in res} if key == "gamma" else lh.to_gamma_numbering(res, res_g)
    elem = {r[0]: lh.element_of(gmap.get(r[0])) for r in res}
    xyz = np.array([a[3] for a in atoms])
    a_res = np.array([a[0] for a in atoms])
    a_name = {a[0]: a[1] for a in atoms}
    helix_names = [e[0] for e in lh.ELEMENTS if e[3] == "helix"]
    rows = []
    for _, c in cons.iterrows():
        d = np.linalg.norm(xyz - np.array([c.cx, c.cy, c.cz]), axis=1)
        near = {}
        for i in np.argsort(d):
            if d[i] > CONTACT:
                break
            near.setdefault(int(a_res[i]), float(d[i]))
        by = {}
        for r, dist in near.items():
            by.setdefault(elem[r], []).append((dist, r))
        es = sorted(by.items(), key=lambda kv: min(kv[1])[0])
        best = (None, np.inf)
        for h in helix_names:
            idx = np.isin(a_res, [r for r in elem if elem[r] == h])
            if idx.any() and d[idx].min() < best[1]:
                best = (h, float(d[idx].min()))
        rows.append({"C": c.C, "pct": round(c.pct, 1), "pct_nr": c.pct_nr, "class": m1.classify(c.pct),
                     "location": " / ".join(e for e, _ in es) or f"no residue within {CONTACT:g} Å",
                     "nearest_helix": best[0], "dist_nearest_helix_A": round(best[1], 2),
                     f"residues_within_{CONTACT:g}A": "; ".join(
                         f"{e} (" + ", ".join(f"{lh.SHORT[a_name[r]]}{r}" for _, r in sorted(v)) + ")"
                         for e, v in es),
                     "reference_based_site": c.nearest_ref_site if str(c.in_reference) == "True" else ""})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- task 2
def bfactors(Wm, cl):
    g = Wm[Wm.label >= 0].groupby("label")["B_norm"].median().rename("median_B_norm")
    t = cl.merge(g, left_on="label", right_index=True)
    t = t[t.n_structures >= MIN_STRUCT].copy()
    t["class"] = t.pct.map(m1.classify)
    rho, p = stats.spearmanr(t.pct, t.median_B_norm)
    hi, ra = t.loc[t["class"] == "high", "median_B_norm"], t.loc[t["class"] == "rare", "median_B_norm"]
    u = stats.mannwhitneyu(hi, ra).pvalue if len(hi) >= 3 and len(ra) >= 3 else np.nan
    out = t[["C", "pct", "n_structures", "median_B_norm", "class"]].round(3)
    return out, {"rho": rho, "p": p, "u_p": u, "med_high": hi.median(), "med_rare": ra.median(),
                 "n_clusters": len(t)}


# ---------------------------------------------------------------- task 3
def resolution(Wm, cl, n_total):
    cons = cl[cl.pct >= LOW]
    cons_labels = set(cons.label)
    res_of = Wm.drop_duplicates("pdb_id").set_index("pdb_id")["resolution_A"]
    n_w = Wm.groupby("pdb_id").size()
    hit = Wm[Wm.label.isin(cons_labels)].groupby("pdb_id")["label"].nunique()
    ps = pd.DataFrame({"resolution_A": res_of, "n_waters_shell": n_w,
                       "conserved_clusters_occupied": hit}).fillna({"conserved_clusters_occupied": 0})
    ps["pct_conserved_occupied"] = 100 * ps.conserved_clusters_occupied / len(cons)
    ps = ps.reset_index().rename(columns={"index": "pdb_id"})
    ok = ps.resolution_A.notna()
    rc, pc = stats.spearmanr(ps.resolution_A[ok], ps.pct_conserved_occupied[ok])
    ra, pa = stats.spearmanr(ps.resolution_A[ok], ps.n_waters_shell[ok])
    # per cluster by resolution group
    grp = pd.cut(res_of, RES_BINS, labels=RES_LABELS)
    members = Wm[Wm.label.isin(cons_labels)].groupby("label")["pdb_id"].apply(set)
    rows = []
    for _, c in cons.iterrows():
        row = {"C": c.C, "pct_all": round(c.pct, 1)}
        for lab in RES_LABELS:
            ids = set(grp[grp == lab].index)
            row[f"pct_{lab}"] = round(100 * len(ids & members[c.label]) / len(ids), 1) if ids else np.nan
            row[f"n_{lab}"] = len(ids)
        rows.append(row)
    no_cons = int((ps.conserved_clusters_occupied == 0).sum())
    return ps, pd.DataFrame(rows), {"rho_res": rc, "p_res": pc, "rho_all": ra, "p_all": pa,
                                    "structures_without_conserved_water": no_cons}


# ---------------------------------------------------------------- figures
def style(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    if grid == "both":
        ax.grid(True, color=GRID, linewidth=0.8)
    elif grid:
        getattr(ax, f"{grid}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def figures(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    keys = ["alpha", "gamma", "delta"]
    order = [e[0] for e in lh.ELEMENTS if e[3] == "helix"]
    rng = np.random.default_rng(0)

    # task 1
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), dpi=300, sharey=True)
    for ax, k in zip(axes, keys):
        t = out[k]["helix"]
        hi = Counter(t[t["class"] == "high"].nearest_helix)
        mo = Counter(t[t["class"] == "moderate"].nearest_helix)
        y = np.arange(len(order))
        ax.barh(y, [mo[h] for h in order], color="#9ec5f4", label="moderate (30–70%)")
        ax.barh(y, [hi[h] for h in order], left=[mo[h] for h in order], color="#256abf",
                label="high (>70%)")
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=9)
        ax.invert_yaxis()
        ax.set_title(SUB[k][0], loc="left", fontsize=11, color=INK)
        ax.set_xlabel("Conserved clusters", fontsize=9, color=INK)
        style(ax, "x")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], fontsize=8, frameon=False, ncol=2, loc="upper right",
               bbox_to_anchor=(0.99, 0.99))
    fig.suptitle("Nearest helix of each conserved Method 1 cluster (≥30% of structures)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(FIG, "m1_sites_per_helix.png"), facecolor="white")
    plt.close(fig)

    # task 2
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), dpi=300, sharey=True)
    for ax, k in zip(axes, keys):
        t, st = out[k]["bf"], out[k]["bf_stats"]
        col = np.where(t["class"] == "high", ACC, MUTED)
        ax.scatter(t.pct, t.median_B_norm, s=9, c=col, alpha=0.7, edgecolors="none", zorder=3)
        ax.axhline(1.0, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.axvline(HIGH, color=GRID, linewidth=1)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 3)
        ax.set_title(f"{SUB[k][0]}   ρ = {st['rho']:.2f}, {p_text(st['p'])}", loc="left",
                     fontsize=10, color=INK)
        ax.set_xlabel("Cluster conservation (% of structures)", fontsize=9, color=INK)
        style(ax)
    axes[0].set_ylabel("Median water B / mean Cα B", fontsize=9, color=INK)
    axes[0].text(2, 1.03, "as mobile as the average Cα", fontsize=7, color=INK2)
    fig.suptitle(f"More conserved clusters hold their water more firmly "
                 f"(each dot = cluster seen in ≥{MIN_STRUCT} structures; blue = >70%)",
                 x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(FIG, "m1_bfactor_vs_conservation.png"), facecolor="white")
    plt.close(fig)

    # task 3 scatter
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2), dpi=300, sharex=True)
    for j, k in enumerate(keys):
        ps, st = out[k]["res"], out[k]["res_stats"]
        x = ps.resolution_A + rng.uniform(-0.02, 0.02, len(ps))
        axes[0, j].scatter(x, ps.pct_conserved_occupied, s=12, color=ACC, alpha=0.55,
                           edgecolors="white", linewidths=0.3, zorder=3)
        axes[0, j].set_ylim(-3, 103)
        axes[0, j].set_title(f"{SUB[k][0]}   ρ = {st['rho_res']:.2f}, {p_text(st['p_res'])}",
                             loc="left", fontsize=10, color=INK)
        axes[1, j].scatter(x, ps.n_waters_shell, s=12, color=MUTED, alpha=0.6,
                           edgecolors="white", linewidths=0.3, zorder=3)
        axes[1, j].set_title(f"{SUB[k][0]}   ρ = {st['rho_all']:.2f}, {p_text(st['p_all'])}",
                             loc="left", fontsize=10, color=INK)
        axes[1, j].set_xlabel("Resolution (Å, lower = better)", fontsize=9, color=INK)
        for i in range(2):
            v = ps.pct_conserved_occupied if i == 0 else ps.n_waters_shell
            step = 0.25
            edges = np.arange(np.floor(ps.resolution_A.min() / step) * step,
                              ps.resolution_A.max() + step, step)
            for lo, hi in zip(edges[:-1], edges[1:]):
                m = (ps.resolution_A >= lo) & (ps.resolution_A < hi)
                if m.sum() >= 3:
                    q1, md, q3 = np.percentile(v[m], [25, 50, 75])
                    mid = (lo + hi) / 2
                    axes[i, j].plot([mid, mid], [q1, q3], color=INK, linewidth=2, zorder=4)
                    axes[i, j].scatter([mid], [md], s=26, color=INK, zorder=5, edgecolors="white")
            style(axes[i, j], "both")
    axes[0, 0].set_ylabel("Conserved clusters (≥30%) with a water\n(% of the subtype's clusters)",
                          fontsize=9, color=INK)
    axes[1, 0].set_ylabel("Waters within 5 Å of the LBD", fontsize=9, color=INK)
    fig.suptitle("Resolution vs conserved clusters occupied (top) and all LBD waters (bottom); "
                 "black = median and middle 50% per 0.25 Å", x=0.01, ha="left", fontsize=10.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(FIG, "m1_resolution_vs_conservation.png"), facecolor="white")
    plt.close(fig)

    # task 3 heatmap
    cmap = LinearSegmentedColormap.from_list(
        "blue", ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
    fig = plt.figure(figsize=(11, 0.27 * HEAT_N + 2.0), dpi=300)
    gs = fig.add_gridspec(1, 3, wspace=0.5, left=0.06, right=0.9, top=0.86, bottom=0.03)
    for j, k in enumerate(keys):
        t = out[k]["res_sites"].head(HEAT_N)
        hx = out[k]["helix"].set_index("C").nearest_helix
        ax = fig.add_subplot(gs[0, j])
        cols = ["pct_all"] + [f"pct_{g}" for g in RES_LABELS]
        M = t[cols].to_numpy(float)
        im = ax.imshow(M, cmap=cmap, vmin=0, vmax=100, aspect="auto")
        for i in range(M.shape[0]):
            for c in range(M.shape[1]):
                if not np.isnan(M[i, c]):
                    ax.text(c, i, f"{M[i, c]:.0f}", ha="center", va="center", fontsize=6.5,
                            color="white" if M[i, c] >= 55 else INK)
        ax.axvline(0.5, color="white", linewidth=4)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["all"] + [f"{g}\nn={t[f'n_{g}'].iloc[0]}" for g in RES_LABELS],
                           fontsize=7.5, color=INK)
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels([f"{c}  {hx.get(c, '')}" for c in t.C], fontsize=7, color=INK)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(SUB[k][0], loc="left", fontsize=11, color=INK, pad=24)
    cax = fig.add_axes([0.925, 0.4, 0.012, 0.4])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("% of structures with a water in the cluster", fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=INK2)
    cb.outline.set_visible(False)
    fig.suptitle(f"Top {HEAT_N} Method 1 clusters by resolution group (label = nearest helix)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.savefig(os.path.join(FIG, "m1_resolution_heatmap.png"), facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(RES, exist_ok=True)
    _, res_g = lh.read_chain(os.path.join(BASE_DIR, SUB["gamma"][2], "1prg.cif"), "A")
    out, stat_rows = {}, []
    for key, (label, ref, raw) in SUB.items():
        Wm, cl, n_total = membership(key, ref)
        cons = cl[cl.C != ""]
        cons = cons[cons.pct >= LOW]
        hx = helices(key, raw, ref, cons, res_g)
        bf, bfs = bfactors(Wm, cl)
        ps, rsites, rs = resolution(Wm, cl, n_total)
        rsites = rsites.sort_values("pct_all", ascending=False)
        out[key] = {"helix": hx, "bf": bf, "bf_stats": bfs, "res": ps, "res_stats": rs,
                    "res_sites": rsites}
        hx.to_csv(os.path.join(RES, f"followup_helices_{key}.csv"), index=False)
        bf.to_csv(os.path.join(RES, f"followup_bfactor_{key}.csv"), index=False)
        ps.to_csv(os.path.join(RES, f"followup_resolution_{key}.csv"), index=False)
        rsites.to_csv(os.path.join(RES, f"followup_resolution_sites_{key}.csv"), index=False)
        hi = hx[hx["class"] == "high"]
        stat_rows.append({
            "subtype": label, "structures": n_total, "conserved_clusters_ge30": len(hx),
            "high_clusters": len(hi),
            "high_nearest_helix": "; ".join(f"{h} {n}" for h, n in Counter(hi.nearest_helix).most_common()),
            "spearman_conservation_vs_Bnorm": round(bfs["rho"], 3), "p_Bnorm": bfs["p"],
            "median_Bnorm_high": round(bfs["med_high"], 3), "median_Bnorm_rare": round(bfs["med_rare"], 3),
            "mannwhitney_p_high_vs_rare": bfs["u_p"],
            "spearman_resolution_vs_conserved_occupied": round(rs["rho_res"], 3), "p_resolution": rs["p_res"],
            "spearman_resolution_vs_all_waters": round(rs["rho_all"], 3),
            "high_clusters_still_gt70_at_le2A": int((rsites.loc[rsites.pct_all > HIGH, "pct_≤2.0 Å"] > HIGH).sum()),
            "structures_without_conserved_water": rs["structures_without_conserved_water"],
        })
        print(f"== {label}: {len(hx)} conserved clusters ({len(hi)} high); "
              f"B_norm ρ = {bfs['rho']:.2f} ({p_text(bfs['p'])}), high {bfs['med_high']:.2f} vs rare "
              f"{bfs['med_rare']:.2f}; resolution ρ = {rs['rho_res']:.2f} ({p_text(rs['p_res'])}), "
              f"all waters ρ = {rs['rho_all']:.2f}")
        print(hx.head(8)[["C", "pct", "location", "nearest_helix", "dist_nearest_helix_A",
                          "reference_based_site"]].to_string(index=False))
    st = pd.DataFrame(stat_rows)
    st.to_csv(os.path.join(RES, "followup_stats.csv"), index=False)
    print(st.drop(columns=["p_Bnorm", "p_resolution", "mannwhitney_p_high_vs_rare"]).to_string(index=False))
    figures(out)
    print("Tables in", RES, "| figures in", FIG)


if __name__ == "__main__":
    main()
