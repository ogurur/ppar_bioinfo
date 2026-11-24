# =====================================================================
# Method 1, step 2: cluster ALL waters, conservation per cluster, sensitivity
# ---------------------------------------------------------------------
# Waters of every structure are placed in the main reference frame with the
# core-helix superposition from method1_prepare.py, then clustered:
#   greedy  - distance-cutoff clustering (section 4.4 A), seeded from the
#             densest water so a cluster centres on a site, re-centred once
#   dbscan  - DBSCAN (section 4.4 B), min_samples = 3
# Conservation of a cluster = structures with >= 1 water in it / structures
# analysed (section 4.5).  Structures analysed = core RMSD <= MAX_RMSD and at
# least one LBD water (structures with no modelled waters carry no information).
#
# Sensitivity (section 6.2): clustering cutoff 1.0-3.5 Å for both methods;
# resolution filter (all / ≤2.5 / ≤2.0 Å); alignment scope (core vs full LBD);
# including the RMSD outliers; 3 references per subtype.  Sets of high-
# conservation clusters are compared by Jaccard similarity (two clusters are
# the same site when their centroids are within MATCH Å).
#
# Input:  data/processed/method1/*  (method1_prepare.py)
#         data/results/water_sites_conservation_<subtype>.csv (reference-based sites)
# Output: data/results/method1/clusters_<subtype>.csv       the conservation atlas
#         data/results/method1/sensitivity_*.csv
#         figures/method1_*.png, figures/method1_clusters_<subtype>.pdb/.pml
#
# Run:  python method1_cluster.py            (all subtypes; or: python method1_cluster.py . gamma)
#       large subtypes can be split:  python method1_cluster.py . gamma sweep  then  ... gamma sens
#       Needs: numpy, pandas, scipy, scikit-learn, matplotlib
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
from sklearn.cluster import DBSCAN

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(BASE_DIR, "data", "processed", "method1")
RES = os.path.join(BASE_DIR, "data", "results", "method1")
FIG = os.path.join(BASE_DIR, "figures")

MAIN = ("greedy", 1.5)                     # main setting: method, cutoff (Å)
CUTOFFS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]   # manual: 1.5-3.5; 1.0 added for context
MAX_RMSD = 2.0                             # core Cα RMSD above this = outlier
MIN_SAMPLES = 3                            # DBSCAN
MATCH = 1.5                                # Å, same site when centroids this close
HIGH, LOW = 70.0, 30.0
CELL_TOL = 0.02                            # near-duplicates: same lab + space group, cell within 2%

SUBTYPES = {"alpha": ("PPARα", "1k7l"), "gamma": ("PPARγ", "1prg"), "delta": ("PPARδ", "1gwx")}


# ---------------------------------------------------------------------
def load(key):
    W = pd.read_csv(os.path.join(PRE, f"waters_{key}.csv"), dtype={"water_resi": str})
    T = pd.read_csv(os.path.join(PRE, f"transforms_{key}.csv"))
    return W, T


def place(W, T, reference, scope):
    """Water coordinates moved onto `reference` with the given superposition."""
    t = T[(T.reference == reference) & (T.scope == scope)].set_index("pdb_id")
    out, keep = np.zeros((len(W), 3)), np.zeros(len(W), bool)
    for pid, idx in W.groupby("pdb_id").indices.items():
        if pid not in t.index or pd.isna(t.loc[pid, "R00"]):
            continue
        r = t.loc[pid]
        R = np.array([[r[f"R{i}{j}"] for j in range(3)] for i in range(3)])
        tt = np.array([r.t0, r.t1, r.t2])
        out[idx] = W.iloc[idx][["x", "y", "z"]].to_numpy() @ R.T + tt
        keep[idx] = True
    return out, keep, t["rmsd"]


def greedy(X, cutoff):
    """Density-seeded distance-cutoff clustering -> labels."""
    tree = cKDTree(X)
    dens = np.array([len(n) for n in tree.query_ball_point(X, cutoff)])
    labels = -np.ones(len(X), int)
    cid = 0
    for i in np.argsort(-dens, kind="stable"):
        if labels[i] != -1:
            continue
        mem = [j for j in tree.query_ball_point(X[i], cutoff) if labels[j] == -1]
        cen = X[mem].mean(0)                                   # re-centre once
        mem = [j for j in tree.query_ball_point(cen, cutoff) if labels[j] == -1]
        if i not in mem:
            mem.append(i)
        labels[mem] = cid
        cid += 1
    return labels


def cluster(X, method, cutoff):
    if method == "greedy":
        return greedy(X, cutoff)
    return DBSCAN(eps=cutoff, min_samples=MIN_SAMPLES).fit(X).labels_


def summarise(labels, X, pids, n_total, extra=None, wmap=None):
    """One row per cluster (noise excluded)."""
    df = pd.DataFrame({"label": labels, "pdb_id": pids, "x": X[:, 0], "y": X[:, 1], "z": X[:, 2]})
    if extra is not None:
        for k, v in extra.items():
            df[k] = v
    df = df[df.label >= 0]
    rows = []
    for lab, g in df.groupby("label"):
        c = g[["x", "y", "z"]].to_numpy().mean(0)
        d = np.linalg.norm(g[["x", "y", "z"]].to_numpy() - c, axis=1)
        per_struct = g.groupby("pdb_id").size()
        row = {"label": lab, "n_structures": per_struct.size, "n_waters": len(g),
               "pct": 100 * per_struct.size / n_total, "cx": c[0], "cy": c[1], "cz": c[2],
               "radius90": float(np.percentile(d, 90)),
               "pct_structures_multi": 100 * (per_struct > 1).mean()}
        if wmap is not None:
            row["pct_nr"] = 100 * sum(wmap.get(p, 0) for p in per_struct.index) / wmap["__total__"]
        if extra is not None and "B_norm" in g:
            row["median_B_norm"] = float(g["B_norm"].median())
            row["resolutions"] = g.drop_duplicates("pdb_id")["resolution_A"].to_numpy()
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values("pct", ascending=False).reset_index(drop=True) if len(out) else out


def high_set(cl, thr=HIGH):
    return cl.loc[cl.pct > thr, ["cx", "cy", "cz"]].to_numpy() if len(cl) else np.zeros((0, 3))


def match_sets(A, B, d=MATCH):
    """Greedy one-to-one matching of two centroid sets -> number of matches."""
    if not len(A) or not len(B):
        return 0
    D = np.linalg.norm(A[:, None] - B[None], axis=2)
    n = 0
    while True:
        i, j = np.unravel_index(np.argmin(D), D.shape)
        if D[i, j] > d:
            return n
        n += 1
        D[i, :] = np.inf
        D[:, j] = np.inf


def jaccard(A, B):
    m = match_sets(A, B)
    u = len(A) + len(B) - m
    return m / u if u else np.nan


def classify(p):
    return "high" if p > HIGH else ("moderate" if p >= LOW else "rare")


def near_duplicate_groups(pdb_ids):
    """Section 3.3: group structures from the same depositing group, same space
    group and unit-cell lengths within CELL_TOL.  -> {pdb_id: group id}"""
    meta = pd.read_csv(os.path.join(BASE_DIR, "data", "metadata.csv"))
    meta["pdb_id"] = meta["pdb_id"].str.lower()
    m = meta.set_index("pdb_id").reindex(sorted(pdb_ids))
    parent = {p: p for p in m.index}

    def find(p):
        while parent[p] != p:
            parent[p] = parent[parent[p]]
            p = parent[p]
        return p

    for _, g in m.groupby([m.depositing_group.fillna("?"), m.space_group.fillna("?")]):
        ids = list(g.index)
        cell = g[["cell_a", "cell_b", "cell_c"]].to_numpy(float)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                if np.all(np.abs(cell[i] - cell[j]) <= CELL_TOL * np.maximum(cell[i], cell[j])):
                    parent[find(ids[i])] = find(ids[j])
    roots = {p: find(p) for p in m.index}
    names = {r: f"G{i + 1}" for i, r in enumerate(sorted(set(roots.values())))}
    return {p: names[r] for p, r in roots.items()}


# ---------------------------------------------------------------------
def run_subtype(key, label, ref, refs_all, stage="all"):
    W, T = load(key)
    X, ok, rmsd = place(W, T, ref, "core")
    good = W.pdb_id.map(lambda p: rmsd.get(p, np.inf) <= MAX_RMSD).to_numpy() & ok
    Wm, Xm = W[good].reset_index(drop=True), X[good]
    n_total = Wm.pdb_id.nunique()
    n_out = int((rmsd > MAX_RMSD).sum())
    print(f"== {label}: {n_total} structures with waters and core RMSD ≤ {MAX_RMSD} Å "
          f"({n_out} outliers left out), {len(Wm)} waters")

    extra = {"B_norm": Wm.B_norm.to_numpy(), "resolution_A": Wm.resolution_A.to_numpy()}
    grp = near_duplicate_groups(Wm.pdb_id.unique())
    gsize = pd.Series(grp).value_counts()
    wmap = {p: 1.0 / gsize[g] for p, g in grp.items()}
    wmap["__total__"] = float(len(gsize))            # = number of independent groups
    pd.DataFrame({"pdb_id": list(grp), "group": list(grp.values())}) \
        .assign(group_size=lambda d: d.group.map(gsize)) \
        .sort_values(["group_size", "group"], ascending=[False, True]) \
        .to_csv(os.path.join(RES, f"near_duplicate_groups_{key}.csv"), index=False)
    print(f"   near-duplicate groups: {len(gsize)} independent groups for {n_total} structures; "
          f"largest {gsize.max()} ({gsize.index[0]})")

    # --- main clustering + cutoff sweep --------------------------------
    sweep, cls = [], {}
    for method in (("greedy", "dbscan") if stage in ("all", "sweep") else ()):
        for cut in CUTOFFS:
            cl = summarise(cluster(Xm, method, cut), Xm, Wm.pdb_id.to_numpy(), n_total, extra)
            cls[(method, cut)] = cl
    main = summarise(cluster(Xm, *MAIN), Xm, Wm.pdb_id.to_numpy(), n_total, extra, wmap)
    cls[MAIN] = main
    if stage == "sens":
        cls = {}
    for (method, cut), cl in cls.items():
        hs = high_set(cl)
        sweep.append({"subtype": label, "method": method, "cutoff_A": cut,
                      "n_clusters_ge2_structures": int((cl.n_structures >= 2).sum()),
                      "n_high": int((cl.pct > HIGH).sum()),
                      "n_moderate": int(((cl.pct >= LOW) & (cl.pct <= HIGH)).sum()),
                      "median_radius90_high": float(cl.loc[cl.pct > HIGH, "radius90"].median())
                      if (cl.pct > HIGH).any() else np.nan,
                      "median_pct_multi_high": float(cl.loc[cl.pct > HIGH, "pct_structures_multi"].median())
                      if (cl.pct > HIGH).any() else np.nan,
                      "jaccard_high_vs_main": round(jaccard(hs, high_set(main)), 3)})

    # --- atlas: main clusters (>= LOW %) with reference-based site match -----
    sites = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                     f"water_sites_conservation_{key}.csv"))
    S = sites[["ref_x", "ref_y", "ref_z"]].to_numpy()
    at = main[(main.pct >= LOW) | (main.pct_nr >= LOW)].copy()
    at.insert(0, "C", [f"C{i + 1}" for i in range(len(at))])
    d = np.linalg.norm(at[["cx", "cy", "cz"]].to_numpy()[:, None] - S[None], axis=2)
    at["nearest_ref_site"] = sites["W"].to_numpy()[d.argmin(1)]
    at["nearest_ref_site_pct"] = sites["pct_1.0A"].to_numpy()[d.argmin(1)]
    at["dist_to_ref_site_A"] = d.min(1).round(2)
    at["in_reference"] = at["dist_to_ref_site_A"] <= MATCH
    at["class"] = at.pct.map(classify)
    at["class_nr"] = at.pct_nr.map(classify)
    # conservation of the same cluster among better-resolved structures only
    res_of = Wm.drop_duplicates("pdb_id").set_index("pdb_id")["resolution_A"]
    for cutres in (2.5, 2.0):
        n_sub = int((res_of <= cutres).sum())
        at[f"pct_res_le_{cutres}"] = [100 * (r <= cutres).sum() / n_sub for r in at["resolutions"]]
    at = at.drop(columns=["label", "resolutions"])
    for c in ("pct", "pct_nr", "cx", "cy", "cz", "radius90", "pct_structures_multi", "median_B_norm",
              "pct_res_le_2.5", "pct_res_le_2.0"):
        at[c] = at[c].round(2)

    # --- sensitivity: resolution, scope, outliers, references -----------------
    sens = []
    hs_main = high_set(main)

    def add(test, setting, cl, n):
        sens.append({"subtype": label, "test": test, "setting": setting, "n_structures": n,
                     "n_high": int((cl.pct > HIGH).sum()) if len(cl) else 0,
                     "n_moderate": int(((cl.pct >= LOW) & (cl.pct <= HIGH)).sum()) if len(cl) else 0,
                     "jaccard_high_vs_main": round(jaccard(high_set(cl), hs_main), 3)})

    if stage == "sweep":
        refs_all = []
    add("main", f"{MAIN[0]} {MAIN[1]} Å, core, {ref.upper()}, all resolutions", main, n_total)
    nr = main.assign(pct=main.pct_nr)
    add("near-duplicates", f"each lab/crystal-form group counted once ({len(gsize)} groups)",
        nr, len(gsize))
    for cutres in ((2.5, 2.0) if stage != "sweep" else ()):
        m = (Wm.resolution_A <= cutres).to_numpy()
        n = Wm[m].pdb_id.nunique()
        cl = summarise(cluster(Xm[m], *MAIN), Xm[m], Wm[m].pdb_id.to_numpy(), n)
        add("resolution", f"≤{cutres} Å", cl, n)
    Xf, okf, rf = place(W, T, ref, "full")
    if stage != "sweep":
        gf = W.pdb_id.map(lambda p: rf.get(p, np.inf) <= 99).to_numpy() & okf & good
        cl = summarise(cluster(Xf[gf], *MAIN), Xf[gf], W[gf].pdb_id.to_numpy(), W[gf].pdb_id.nunique())
        add("alignment scope", "full LBD Cα", cl, W[gf].pdb_id.nunique())
        ga = ok.copy()
        cl = summarise(cluster(X[ga], *MAIN), X[ga], W[ga].pdb_id.to_numpy(), W[ga].pdb_id.nunique())
        add("outliers", f"including RMSD > {MAX_RMSD} Å", cl, W[ga].pdb_id.nunique())

    # alternative references: cluster in their frame, map centroids back
    tm = T[(T.reference == ref) & (T.scope == "core")].set_index("pdb_id")
    main_mod = main[main.pct >= LOW]
    for alt in refs_all:
        if alt == ref:
            continue
        Xa, oka, ra = place(W, T, alt, "core")
        ga = oka & W.pdb_id.map(lambda p: ra.get(p, np.inf) <= MAX_RMSD).to_numpy()
        n = W[ga].pdb_id.nunique()
        cl = summarise(cluster(Xa[ga], *MAIN), Xa[ga], W[ga].pdb_id.to_numpy(), n)
        r = tm.loc[alt]
        R = np.array([[r[f"R{i}{j}"] for j in range(3)] for i in range(3)])
        tt = np.array([r.t0, r.t1, r.t2])
        C = cl[["cx", "cy", "cz"]].to_numpy() @ R.T + tt
        cl[["cx", "cy", "cz"]] = C
        add("reference", alt.upper(), cl, n)
        # conservation change for the main sites >= LOW %
        D = np.linalg.norm(main_mod[["cx", "cy", "cz"]].to_numpy()[:, None] - C[None], axis=2)
        j = D.argmin(1)
        hit = D.min(1) <= MATCH
        diff = np.abs(main_mod.pct.to_numpy()[hit] - cl.pct.to_numpy()[j][hit])
        rho = stats.spearmanr(main_mod.pct.to_numpy()[hit], cl.pct.to_numpy()[j][hit])[0] \
            if hit.sum() > 3 else np.nan
        sens[-1].update({"main_sites_ge30_found": f"{hit.sum()}/{len(main_mod)}",
                         "median_abs_change_pct": round(float(np.median(diff)), 1) if len(diff) else np.nan,
                         "max_abs_change_pct": round(float(diff.max()), 1) if len(diff) else np.nan,
                         "sites_changed_gt10pct": int((diff > 10).sum()),
                         "spearman_rank_vs_main": round(rho, 3)})

    rms = pd.DataFrame({"pdb_id": rmsd.index, "rmsd_core": rmsd.values})
    rms["rmsd_full"] = rms.pdb_id.map(rf)
    rms["outlier"] = rms.rmsd_core > MAX_RMSD
    outs = [("clusters", at), ("rmsd", rms)]
    if stage in ("all", "sweep"):
        outs.append(("sweep", pd.DataFrame(sweep)))
    if stage in ("all", "sens"):
        outs.append(("sens", pd.DataFrame(sens)))
    for k, v in outs:
        v.to_csv(os.path.join(RES, f"{k}_{key}.csv"), index=False)
    write_pdb(key, ref, at)
    print(f"   main ({MAIN[0]} {MAIN[1]} Å): {int((at.pct > HIGH).sum())} high, "
          f"{int(((at.pct >= LOW) & (at.pct <= HIGH)).sum())} moderate clusters; "
          f"{int((~at.in_reference).sum())} of them not at a reference water")
    print(f"   non-redundant: {int((at.pct_nr > HIGH).sum())} high, "
          f"{int(((at.pct_nr >= LOW) & (at.pct_nr <= HIGH)).sum())} moderate")
    print(at.head(10)[["C", "pct", "pct_nr", "n_structures", "radius90", "pct_structures_multi",
                       "median_B_norm", "nearest_ref_site", "dist_to_ref_site_A",
                       "pct_res_le_2.0"]].to_string(index=False))
    return


def write_pdb(key, ref, at):
    """Cluster centroids as pseudo-waters, B-factor column = conservation %."""
    lines = []
    for i, r in enumerate(at.itertuples(), 1):
        lines.append("HETATM%5d  O   HOH W%4d    %8.3f%8.3f%8.3f  1.00%6.2f           O"
                     % (i, i, r.cx, r.cy, r.cz, r.pct))
    with open(os.path.join(FIG, f"method1_clusters_{key}.pdb"), "w") as f:
        f.write("\n".join(lines) + "\nEND\n")
    folder = {"alpha": "PPARalpha", "gamma": "PPARgamma", "delta": "PPARdelta"}[key]
    pml = f"""# Method 1 clusters (>=30%) on the {ref.upper()} reference; sphere colour/size = conservation
# run from the project folder:  @figures/method1_clusters_{key}.pml
reinitialize
load data/raw/{folder}/{ref}.cif, ref
remove ref and not chain A
remove ref and resn HOH
hide everything
show cartoon, ref
color grey80, ref
bg_color white
load figures/method1_clusters_{key}.pdb, sites
show spheres, sites
set sphere_scale, 0.35, sites
set sphere_scale, 0.55, sites and b > 70
spectrum b, yellow_red, sites, minimum=30, maximum=90
label sites and b > 70, "C%s" % resi.strip()
set label_color, black
orient ref
"""
    with open(os.path.join(FIG, f"method1_clusters_{key}.pml"), "w") as f:
        f.write(pml)


# ---------------------------------------------------------------------
INK, INK2, GRID, ACC, MUTED = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#a8a69f"


def style(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    if grid:
        getattr(ax, f"{grid}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plots():
    figs = {}
    for k, (label, _) in SUBTYPES.items():
        figs[k] = (label, pd.read_csv(os.path.join(RES, f"rmsd_{k}.csv")),
                   pd.read_csv(os.path.join(RES, f"sweep_{k}.csv")),
                   pd.read_csv(os.path.join(RES, f"clusters_{k}.csv")))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = ["alpha", "gamma", "delta"]

    # 1 - RMSD core vs full
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), dpi=300, sharey=False)
    bins = np.arange(0, 8.01, 0.2)
    for ax, k in zip(axes, keys):
        label, rms, _, _ = figs[k]
        ax.hist(rms.rmsd_full.dropna(), bins=bins, color=MUTED, label="full LBD Cα")
        ax.hist(rms.rmsd_core.dropna(), bins=bins, color=ACC, alpha=0.85, label="core-helix Cα")
        ax.axvline(MAX_RMSD, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.set_title(f"{label}: core median {rms.rmsd_core.median():.2f} Å, "
                     f"{int(rms.outlier.sum())} > {MAX_RMSD:g} Å", loc="left", fontsize=9.5, color=INK)
        ax.set_xlabel("RMSD to reference (Å)", fontsize=9, color=INK)
        style(ax)
    axes[0].set_ylabel("Structures", fontsize=9, color=INK)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Core-helix superposition fits much more tightly than the full LBD",
                 x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(FIG, "method1_rmsd.png"), facecolor="white")
    plt.close(fig)

    # 2 - cutoff sensitivity
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), dpi=300, sharex=True)
    for j, k in enumerate(keys):
        label, _, sw, _ = figs[k]
        for method, col in (("greedy", ACC), ("dbscan", INK2)):
            s = sw[sw.method == method]
            axes[0, j].plot(s.cutoff_A, s.n_high, marker="o", color=col, linewidth=2, markersize=4,
                            label=f"{method} (high)")
            axes[0, j].plot(s.cutoff_A, s.n_moderate, marker="o", color=col, linewidth=1.2,
                            markersize=3, linestyle=(0, (3, 2)), label=f"{method} (moderate)")
            axes[1, j].plot(s.cutoff_A, s.jaccard_high_vs_main, marker="o", color=col, linewidth=2,
                            markersize=4, label=method)
        axes[0, j].set_title(label, loc="left", fontsize=11, color=INK)
        axes[1, j].set_ylim(-0.03, 1.03)
        axes[1, j].set_xlabel("Clustering cutoff (Å)", fontsize=9, color=INK)
        for i in range(2):
            axes[i, j].axvline(MAIN[1], color=GRID, linewidth=6, zorder=0)
            style(axes[i, j])
    axes[0, 0].set_ylabel("Number of clusters", fontsize=9, color=INK)
    axes[1, 0].set_ylabel(f"Jaccard of high sites vs\nmain ({MAIN[0]}, {MAIN[1]:g} Å)", fontsize=9, color=INK)
    axes[0, 2].legend(fontsize=7, frameon=False, loc="upper right")
    fig.suptitle("Clustering sensitivity: high (solid) and moderate (dashed) clusters by cutoff; "
                 "grey band = main setting", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(FIG, "method1_sensitivity_cutoff.png"), facecolor="white")
    plt.close(fig)

    # 3 - ranked clusters, new sites highlighted
    n_show = 30
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.4), dpi=300, sharey=True)
    for ax, k in zip(axes, keys):
        label, _, _, at = figs[k]
        a = at.head(n_show)
        x = np.arange(len(a))
        col = np.where(a.in_reference, MUTED, ACC)
        ax.bar(x, a.pct, color=col, width=0.75)
        ax.scatter(x, a.pct_nr, marker="D", s=12, color=INK, zorder=4,
                   label="non-redundant (lab/crystal-form groups counted once)")
        for y, txt in ((HIGH, "70%"), (LOW, "30%")):
            ax.axhline(y, color=INK2, linewidth=1, linestyle=(0, (4, 3)))
            ax.text(len(a) - 0.4, y + 1.5, txt, ha="right", fontsize=7, color=INK2)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c}\n{w}" if r else c for c, w, r in
                            zip(a.C, a.nearest_ref_site, a.in_reference)], fontsize=6.5, color=INK)
        ax.set_xlim(-0.6, n_show - 0.4)
        ax.set_ylim(0, 100)
        new = int((~at.in_reference).sum())
        ax.set_title(f"{label}: high {int((at.pct > HIGH).sum())} (non-redundant "
                     f"{int((at.pct_nr > HIGH).sum())}); {new} of the clusters ≥30% are "
                     f"not at any reference water (blue)", loc="left", fontsize=9.5, color=INK)
        style(ax)
    axes[1].set_ylabel("Structures with a water in the cluster (%)", fontsize=9, color=INK)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.suptitle(f"Method 1 clusters ({MAIN[0]}, {MAIN[1]:g} Å), ranked; "
                 "label under a grey bar = matching reference-based site", x=0.01, ha="left",
                 fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(FIG, "method1_clusters_ranked.png"), facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    refs = pd.read_csv(os.path.join(PRE, "references.csv"))
    only = sys.argv[2].split(",") if len(sys.argv) > 2 else list(SUBTYPES)
    stage = sys.argv[3] if len(sys.argv) > 3 else "all"     # all | sweep | sens
    for key in only:
        if key in SUBTYPES:
            label, ref = SUBTYPES[key]
            run_subtype(key, label, ref, refs[refs.subtype == key].reference.tolist(), stage)
    done = [k for k in SUBTYPES if all(os.path.isfile(os.path.join(RES, f"{f}_{k}.csv"))
                                       for f in ("sens", "sweep"))]
    if len(done) < len(SUBTYPES):
        print("still to run:", [k for k in SUBTYPES if k not in done])
        return
    pd.concat(pd.read_csv(os.path.join(RES, f"sweep_{k}.csv")) for k in SUBTYPES) \
        .to_csv(os.path.join(RES, "sensitivity_cutoff.csv"), index=False)
    sens = pd.concat(pd.read_csv(os.path.join(RES, f"sens_{k}.csv")) for k in SUBTYPES)
    sens.to_csv(os.path.join(RES, "sensitivity_other.csv"), index=False)
    print(sens.to_string(index=False))
    plots()
    print("Tables in", RES, "| figures in", FIG)


if __name__ == "__main__":
    main()
