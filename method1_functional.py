# =====================================================================
# Are conserved water clusters associated with functional sites?
# (project description, sections 7.2 and 7.3)
# ---------------------------------------------------------------------
# For every Method 1 cluster seen in >= MIN_STRUCT structures:
#   * distance from the cluster centre to functional regions in the reference:
#       pocket     - ligand-binding pocket residues R288, S289, H323, K367, H449, Y473
#       AF-2       - H12 (466-477), charge clamp K301/R304, and PPARγ residues that
#                    touch the coactivator peptide in 1FM6
#       RXR        - PPARγ residues that touch RXRα in 1FM6 (heterodimer interface)
#       selectivity- T253 and I354 as listed in the project description (check numbering)
#     all in PPARγ1 numbering, carried to PPARα / PPARδ by sequence alignment;
#     "near" = within NEAR Å of any heavy atom of the region.
#   * ligand: in ligand-bound structures, how often a ligand atom sits on the
#     site (≤ OCCUPY Å) and how often the cluster's water touches the ligand
#     (≤ TOUCH Å); crystallisation additives and ions are not ligands.
#   * apo vs bound occupancy and the manual's classes (ligand_displaced,
#     ligand_recruited, universally_conserved, state_dependent).
#   * burial: protein heavy atoms within 8 Å of the centre.
# Statistics: Fisher exact test, high (>70 %) vs rare (<30 %) clusters, for
# being near each region; Benjamini-Hochberg correction over all tests.
#
# Input:  data/processed/method1/*, data/results/method1/clusters_<subtype>.csv,
#         data/raw/..., method1_cluster.py, label_helices.py (imported)
# Output: data/results/method1/functional_<subtype>.csv, functional_stats.csv,
#         functional_residues.csv; data/processed/method1/ligands_<subtype>.csv (cache);
#         figures/m1_functional_regions.png, m1_apo_vs_bound.png,
#         figures/m1_functional_<subtype>.pml
# Run:  python method1_functional.py            (or: python method1_functional.py . gamma)
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
ONLY = sys.argv[2].split(",") if len(sys.argv) > 2 else None
sys.path.insert(0, BASE_DIR)
sys.argv = sys.argv[:1] + [BASE_DIR]
import method1_cluster as m1                   # noqa: E402
import label_helices as lh                     # noqa: E402

PRE = os.path.join(BASE_DIR, "data", "processed", "method1")
RES = os.path.join(BASE_DIR, "data", "results", "method1")
FIG = os.path.join(BASE_DIR, "figures")

MIN_STRUCT = 3
NEAR = 5.0          # Å, cluster centre to any heavy atom of a region
OCCUPY = 2.0        # Å, ligand atom on the water site
TOUCH = 3.5         # Å, cluster water to ligand atom
IFACE = 4.5         # Å, residue contact with RXR / coactivator in 1FM6
MIN_LIG_ATOMS = 8
MIN_APO = 10        # apo vs bound only for subtypes with at least this many apo structures
HIGH, LOW = m1.HIGH, m1.LOW

SUB = {"alpha": ("PPARα", "1k7l", "data/raw/PPARalpha"),
       "gamma": ("PPARγ", "1prg", "data/raw/PPARgamma"),
       "delta": ("PPARδ", "1gwx", "data/raw/PPARdelta")}

POCKET_G = {288, 289, 323, 367, 449, 473}
CLAMP_G = {301, 304}
H12_G = set(range(466, 478))
SELECT_G = {253, 354}
COMPLEX = ("data/raw/PPARgamma/1fm6.cif", ["D", "X"])   # PPARγ chains of 1FM6

WATER = {"HOH", "WAT", "DOD"}
POLYMER = set(lh.AA) | {"CME", "CSO", "CSD", "SEP", "TPO", "PTR", "KCX", "OCS", "CSX", "MLY", "ACE", "NH2"}
ADDITIVES = {
    "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "2PE", "12P", "15P", "PE4", "PE5", "PE8",
    "XPE", "7PE", "P4G", "MPD", "MRD", "PGO", "B7G", "BOG", "LDA", "C8E", "LMT", "LMU", "DDQ",
    "CPS", "MES", "EPE", "TRS", "B3P", "BTB", "CIT", "FLC", "TLA", "MLI", "MLA", "ACT", "ACY",
    "DMS", "SO4", "PO4", "CL", "NA", "K", "MG", "ZN", "CA", "IOD", "BR", "NI", "CD", "MN", "CO",
    "FMT", "IPA", "EOH", "MOH", "BME", "DTT", "TCE", "NO3", "NH4", "SCN", "AZI", "CO3", "BCT",
    "HEZ", "IMD", "NHE", "SPM", "PUT", "URE", "GAI", "TAM", "BU3", "BU1", "PEO", "CXS", "SIN",
    "MAL", "SUC", "GLC", "BGC", "TRE", "EDK"}


# ---------------------------------------------------------------- reading
def read_atoms(path):
    """All first-model, first-altloc heavy atoms: DataFrame chain, resi, resn, name, group, x, y, z."""
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
    rows = [r for r in rows if len(r) == len(tags)]
    if "pdbx_PDB_model_num" in c and rows:
        m0 = rows[0][c["pdbx_PDB_model_num"]]
        rows = [r for r in rows if r[c["pdbx_PDB_model_num"]] == m0]
    if "label_alt_id" in c:
        rows = [r for r in rows if r[c["label_alt_id"]] in (".", "?", "A", "1")]
    rows = [r for r in rows if r[c["type_symbol"]] not in ("H", "D")]
    rn = c.get("auth_comp_id", c["label_comp_id"])
    an = c.get("auth_atom_id", c["label_atom_id"])
    return pd.DataFrame({
        "chain": [r[c["auth_asym_id"]] for r in rows], "resi": [r[c["auth_seq_id"]] for r in rows],
        "resn": [r[rn] for r in rows], "name": [r[an] for r in rows],
        "x": [float(r[c["Cartn_x"]]) for r in rows], "y": [float(r[c["Cartn_y"]]) for r in rows],
        "z": [float(r[c["Cartn_z"]]) for r in rows]})


def ligands_of(atoms, chain):
    """Ligand heavy atoms within 5 Å of the chain's protein atoms."""
    prot = atoms[(atoms.chain == chain) & atoms.resn.isin(POLYMER)][["x", "y", "z"]].to_numpy()
    het = atoms[~atoms.resn.isin(POLYMER | WATER | ADDITIVES)]
    if not len(het) or not len(prot):
        return het.iloc[0:0]
    size = het.groupby(["chain", "resn", "resi"]).name.transform("size")
    het = het[size >= MIN_LIG_ATOMS]
    if not len(het):
        return het
    d, _ = cKDTree(prot).query(het[["x", "y", "z"]].to_numpy())
    het = het.assign(d=d)
    keep = het.groupby(["chain", "resn", "resi"]).d.transform("min") <= 5.0
    return het[keep].drop(columns="d")


def ligand_table(key, raw, pdb_chain):
    """Cache of ligand atoms (original coordinates) for every analysed structure."""
    path = os.path.join(PRE, f"ligands_{key}.csv")
    if os.path.isfile(path):
        return pd.read_csv(path, dtype={"resi": str})
    out = []
    for pid, ch in pdb_chain.items():
        f = os.path.join(BASE_DIR, raw, f"{pid}.cif")
        if os.path.isfile(f):
            lig = ligands_of(read_atoms(f), ch)
            if len(lig):
                out.append(lig.assign(pdb_id=pid))
    L = pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["chain", "resi", "resn", "name", "x", "y", "z", "pdb_id"])
    L.to_csv(path, index=False)
    return L


def transform(T, ref, pid, xyz):
    t = T[(T.reference == ref) & (T.scope == "core") & (T.pdb_id == pid)]
    if not len(t) or pd.isna(t.iloc[0]["R00"]):
        return None
    r = t.iloc[0]
    R = np.array([[r[f"R{i}{j}"] for j in range(3)] for i in range(3)])
    return xyz @ R.T + np.array([r.t0, r.t1, r.t2])


# ---------------------------------------------------------------- regions
def interface_residues_gamma(res_g):
    """PPARγ1 residue numbers touching RXR and coactivator chains in 1FM6,
    kept only if found for both PPARγ copies (removes crystal contacts)."""
    path = os.path.join(BASE_DIR, COMPLEX[0])
    atoms = read_atoms(path)
    ch = pd.read_csv(os.path.join(BASE_DIR, "lbd_chains_PPARgamma.csv"), dtype=str).fillna("")
    note = ch.loc[ch.pdb_id == "1fm6", "note"].iloc[0]
    desc = dict(x.split("=", 1) for x in note.replace("skipped: ", "").split(", ") if "=" in x)
    rxr = [c for c, d in desc.items() if "RXR" in d.upper()]
    coa = [c for c, d in desc.items() if "COACTIVATOR" in d.upper() or "NCOA" in d.upper()]
    found = {"RXR": [], "coactivator": []}
    for p in COMPLEX[1]:
        pa = atoms[(atoms.chain == p) & atoms.resn.isin(lh.AA)]
        _, res_p = lh.read_chain(path, p)
        gmap = lh.to_gamma_numbering(res_p, res_g)
        for name, partners in (("RXR", rxr), ("coactivator", coa)):
            pr = atoms[atoms.chain.isin(partners)][["x", "y", "z"]].to_numpy()
            d, _ = cKDTree(pr).query(pa[["x", "y", "z"]].to_numpy())
            hit = set(pa.resi[d <= IFACE].astype(int))
            found[name].append({gmap[r] for r in hit if r in gmap})
    return {k: set.intersection(*v) for k, v in found.items()}


def region_sets(res_g):
    iface = interface_residues_gamma(res_g)
    return {"pocket": POCKET_G, "AF-2": CLAMP_G | H12_G | iface["coactivator"],
            "RXR interface": iface["RXR"], "selectivity": SELECT_G}


# ---------------------------------------------------------------- per subtype
def run(key, regions_g, res_g):
    label, ref, raw = SUB[key]
    W, T = m1.load(key)
    X, ok, rmsd = m1.place(W, T, ref, "core")
    good = W.pdb_id.map(lambda p: rmsd.get(p, np.inf) <= m1.MAX_RMSD).to_numpy() & ok
    Wm, Xm = W[good].reset_index(drop=True), X[good]
    labels = m1.cluster(Xm, *m1.MAIN)
    Wm["label"] = labels
    structs = Wm.drop_duplicates("pdb_id").set_index("pdb_id")["chain"]
    n_total = len(structs)
    cl = m1.summarise(labels, Xm, Wm.pdb_id.to_numpy(), n_total)
    atlas = pd.read_csv(os.path.join(RES, f"clusters_{key}.csv"))
    d = np.linalg.norm(cl[["cx", "cy", "cz"]].to_numpy()[:, None] - atlas[["cx", "cy", "cz"]].to_numpy()[None], axis=2)
    cl["C"] = np.where(d.min(1) < 0.05, atlas["C"].to_numpy()[d.argmin(1)], "")
    cl = cl[cl.n_structures >= MIN_STRUCT].copy()

    # reference atoms and region atoms
    atoms, res = lh.read_chain(os.path.join(BASE_DIR, raw, f"{ref}.cif"), "A")
    gmap = {r[0]: r[0] for r in res} if key == "gamma" else lh.to_gamma_numbering(res, res_g)
    a_xyz = np.array([a[3] for a in atoms])
    a_res = np.array([a[0] for a in atoms])
    region_xyz, region_res = {}, {}
    for name, gset in regions_g.items():
        rs = sorted(r for r in {a[0] for a in atoms} if gmap.get(r) in gset)
        region_res[name] = rs
        region_xyz[name] = a_xyz[np.isin(a_res, rs)]
    burial_tree = cKDTree(a_xyz)

    # ligands, moved into the reference frame
    L = ligand_table(key, raw, structs)
    lig_xyz = {}
    for pid, g in L.groupby("pdb_id"):
        xyz = transform(T, ref, pid, g[["x", "y", "z"]].to_numpy())
        if xyz is not None:
            lig_xyz[pid] = xyz
    bound = [p for p in structs.index if p in lig_xyz]
    apo = [p for p in structs.index if p not in lig_xyz]
    lig_trees = {p: cKDTree(v) for p, v in lig_xyz.items()}

    rows = []
    for _, c in cl.iterrows():
        cen = np.array([c.cx, c.cy, c.cz])
        row = {"C": c.C, "pct": round(c.pct, 1), "n_structures": c.n_structures,
               "class": m1.classify(c.pct),
               "burial_atoms_8A": len(burial_tree.query_ball_point(cen, 8.0))}
        near = []
        for name, xyz in region_xyz.items():
            dist = float(np.linalg.norm(xyz - cen, axis=1).min()) if len(xyz) else np.nan
            row[f"dist_{name}"] = round(dist, 2)
            if dist <= NEAR:
                near.append(name)
        row["near_regions"] = "; ".join(near)
        mem = Wm[Wm.label == c.label]
        has = set(mem.pdb_id)
        row["pct_apo"] = round(100 * len(has & set(apo)) / len(apo), 1) if apo else np.nan
        row["pct_bound"] = round(100 * len(has & set(bound)) / len(bound), 1) if bound else np.nan
        occ = sum(1 for p in bound if lig_trees[p].query(cen)[0] <= OCCUPY)
        row["pct_bound_ligand_on_site"] = round(100 * occ / len(bound), 1) if bound else np.nan
        mi = Wm.index[Wm.label == c.label]
        touch, n_b = 0, 0
        for i in mi:
            p = Wm.at[i, "pdb_id"]
            if p in lig_trees:
                n_b += 1
                touch += lig_trees[p].query(Xm[i])[0] <= TOUCH
        row["pct_waters_touching_ligand"] = round(100 * touch / n_b, 1) if n_b else np.nan
        a, b = row["pct_apo"], row["pct_bound"]
        row["ligand_state"] = (f"n/a (only {len(apo)} apo structures)" if len(apo) < MIN_APO else
                               "ligand_displaced" if a >= 70 and b < 30 else
                               "ligand_recruited" if a < 30 and b >= 70 else
                               "universally_conserved" if a >= 50 and b >= 50 else "state_dependent")
        rows.append(row)
    t = pd.DataFrame(rows).sort_values("pct", ascending=False)
    t["touches_ligand"] = t.pct_waters_touching_ligand >= 30
    t.to_csv(os.path.join(RES, f"functional_{key}.csv"), index=False)
    res_rows = [{"subtype": label, "region": n, "n_residues": len(r),
                 "residues_own_numbering": " ".join(map(str, r))} for n, r in region_res.items()]
    print(f"== {label}: {len(t)} clusters (≥{MIN_STRUCT} structures); {len(bound)} ligand-bound, "
          f"{len(apo)} apo structures")
    cons = t[t.pct >= LOW]
    print(cons.head(12)[["C", "pct", "near_regions", "dist_pocket", "dist_AF-2", "dist_RXR interface",
                         "pct_apo", "pct_bound", "pct_waters_touching_ligand", "ligand_state"]]
          .to_string(index=False))
    write_pml(key, ref, raw, region_res, t)
    return t, res_rows, len(apo), len(bound)


def write_pml(key, ref, raw, region_res, t):
    sel = lambda rs: "+".join(map(str, rs)) or "none"
    cons = t[t.pct >= LOW]
    lines = [f"# Conserved Method 1 clusters (>=30%) and functional regions on {ref.upper()}",
             f"# run from the project folder:  @figures/m1_functional_{key}.pml",
             "reinitialize", f"load {raw}/{ref}.cif, ref", "remove ref and not chain A",
             "remove ref and resn HOH", "bg_color white", "hide everything", "show cartoon, ref",
             "color grey90, ref", "set cartoon_transparency, 0.2",
             f"select pocket, ref and resi {sel(region_res['pocket'])}",
             f"select af2, ref and resi {sel(region_res['AF-2'])}",
             f"select rxr, ref and resi {sel(region_res['RXR interface'])}",
             f"select selectivity, ref and resi {sel(region_res['selectivity'])}",
             "color lightorange, af2", "color palecyan, rxr",
             "show sticks, pocket and not name N+C+O", "color orange, pocket and elem C",
             "show sticks, selectivity and not name N+C+O", "color violet, selectivity and elem C",
             f"load figures/method1_clusters_{key}.pdb, sites", "hide everything, sites",
             "set sphere_scale, 0.4"]
    for _, r in cons.iterrows():
        if not isinstance(r.C, str) or not r.C.startswith("C"):
            continue
        n = int(r.C[1:])
        col = "red" if r["class"] == "high" else "yelloworange"
        lines += [f"show spheres, sites and resi {n}", f"color {col}, sites and resi {n}"]
        if r.near_regions:
            lines.append(f'label sites and resi {n}, "{r.C}"')
    lines += ["set label_color, black", "deselect", "orient ref",
              "# orange sticks = pocket, light orange = AF-2 / coactivator groove,",
              "# pale cyan = RXR interface, violet = selectivity residues,",
              "# red spheres = high (>70%), yellow-orange = moderate (30-70%), labelled if near a region"]
    with open(os.path.join(FIG, f"m1_functional_{key}.pml"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- statistics
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in reversed(list(enumerate(o, 1))):
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


def tests(tables):
    rows = []
    for key, t in tables.items():
        hi, ra = t[t["class"] == "high"], t[t["class"] == "rare"]
        for region in ("pocket", "AF-2", "RXR interface", "selectivity", "touches ligand"):
            f = (lambda d: d.touches_ligand) if region == "touches ligand" else \
                (lambda d, r=region: d.near_regions.str.contains(r, regex=False))
            a, b = int(f(hi).sum()), len(hi) - int(f(hi).sum())
            c, d = int(f(ra).sum()), len(ra) - int(f(ra).sum())
            orr, p = stats.fisher_exact([[a, b], [c, d]]) if len(hi) and len(ra) else (np.nan, np.nan)
            rows.append({"subtype": SUB[key][0], "region": region,
                         "high_near": f"{a}/{len(hi)}", "pct_high_near": round(100 * a / max(len(hi), 1), 1),
                         "rare_near": f"{c}/{len(ra)}", "pct_rare_near": round(100 * c / max(len(ra), 1), 1),
                         "odds_ratio": orr, "fisher_p": p})
    s = pd.DataFrame(rows)
    s["bh_q"] = bh(s.fisher_p.fillna(1.0))
    return s


# ---------------------------------------------------------------- figures
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3de"
CLASS_COL = {"high": "#104281", "moderate": "#2a78d6", "rare": "#86b6ef"}


def figures(tables, counts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [k for k in ("alpha", "gamma", "delta") if k in tables]
    regions = ["pocket", "AF-2", "RXR interface", "selectivity", "touches ligand"]

    fig, axes = plt.subplots(1, len(keys), figsize=(11, 4.4), dpi=300, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, keys):
        t = tables[k]
        y = np.arange(len(regions))
        for j, cls in enumerate(("high", "moderate", "rare")):
            d = t[t["class"] == cls]
            vals = []
            for r in regions:
                m = d.touches_ligand if r == "touches ligand" else d.near_regions.str.contains(r, regex=False)
                vals.append(100 * m.mean() if len(d) else 0)
            ax.barh(y + (j - 1) * 0.26, vals, height=0.24, color=CLASS_COL[cls],
                    label=f"{cls} (n={len(d)})")
        ax.set_yticks(y)
        ax.set_yticklabels(["near pocket", "near AF-2 / coactivator", "near RXR interface",
                            "near selectivity residues", "water touches ligand"], fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of clusters in the class", fontsize=9, color=INK)
        ax.set_title(SUB[k][0], loc="left", fontsize=11, color=INK)
        ax.legend(fontsize=7.5, frameon=False, loc="lower right")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.xaxis.grid(True, color=GRID)
        ax.set_axisbelow(True)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.suptitle(f"Where conserved water clusters sit relative to functional regions "
                 f"(near = centre within {NEAR:g} Å)", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG, "m1_functional_regions.png"), facecolor="white")
    plt.close(fig)

    akeys = [k for k in keys if counts[k][0] >= MIN_APO]
    if not akeys:
        return
    fig, axes = plt.subplots(1, len(akeys), figsize=(5.6 * len(akeys), 4.8), dpi=300, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, akeys):
        t = tables[k]
        t = t[t.pct >= LOW]
        for cls in ("moderate", "high"):
            d = t[t["class"] == cls]
            ax.scatter(d.pct_apo, d.pct_bound, s=18, color=CLASS_COL[cls], alpha=0.85,
                       edgecolors="white", linewidths=0.4, label=cls, zorder=3)
        ax.plot([0, 100], [0, 100], color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.axvspan(70, 100, ymin=0, ymax=0.3, color="#f9ece7", zorder=0)
        ax.axhspan(70, 100, xmin=0, xmax=0.3, color="#e7f1fb", zorder=0)
        ax.text(98, 3, "ligand-displaced", ha="right", fontsize=7, color=INK2)
        ax.text(2, 97, "ligand-recruited", va="top", fontsize=7, color=INK2)
        for _, r in t[t.ligand_state.isin(["ligand_displaced", "ligand_recruited"])].iterrows():
            ax.annotate(r.C, (r.pct_apo, r.pct_bound), xytext=(4, 2), textcoords="offset points",
                        fontsize=7, color=INK)
        n_apo, n_bound = counts[k]
        ax.set_title(f"{SUB[k][0]}  (apo n={n_apo}, bound n={n_bound})", loc="left", fontsize=10, color=INK)
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_xlabel("% of apo structures with the water", fontsize=9, color=INK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(True, color=GRID)
        ax.set_axisbelow(True)
        ax.tick_params(colors=INK2, labelsize=8)
    axes[0].set_ylabel("% of ligand-bound structures\nwith the water", fontsize=9, color=INK)
    axes[0].legend(fontsize=7.5, frameon=False, loc="center right")
    skipped = [f"{SUB[k][0]} {counts[k][0]}" for k in keys if k not in akeys]
    fig.suptitle("Apo vs ligand-bound occupancy of conserved clusters (≥30%)\n"
                 f"dashed = same in both; not shown (apo structures < {MIN_APO}): {', '.join(skipped)}",
                 x=0.01, ha="left", fontsize=10.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(FIG, "m1_apo_vs_bound.png"), facecolor="white")
    plt.close(fig)


def main():
    _, res_g = lh.read_chain(os.path.join(BASE_DIR, SUB["gamma"][2], "1prg.cif"), "A")
    regions_g = region_sets(res_g)
    print("Region definitions (PPARγ1 numbering):")
    for k, v in regions_g.items():
        print(f"   {k}: {len(v)} residues: {' '.join(map(str, sorted(v)))}")
    keys = ONLY or list(SUB)
    for key in keys:
        t, rr, n_apo, n_b = run(key, regions_g, res_g)
        pd.DataFrame(rr).to_csv(os.path.join(RES, f"functional_residues_{key}.csv"), index=False)
        pd.DataFrame([{"n_apo": n_apo, "n_bound": n_b}]).to_csv(
            os.path.join(RES, f"functional_counts_{key}.csv"), index=False)
    if not all(os.path.isfile(os.path.join(RES, f"functional_{k}.csv")) for k in SUB):
        return
    tables = {k: pd.read_csv(os.path.join(RES, f"functional_{k}.csv")).fillna({"near_regions": ""})
              for k in SUB}
    counts = {k: tuple(pd.read_csv(os.path.join(RES, f"functional_counts_{k}.csv")).iloc[0])
              for k in SUB}
    s = tests(tables)
    s.to_csv(os.path.join(RES, "functional_stats.csv"), index=False)
    pd.concat(pd.read_csv(os.path.join(RES, f"functional_residues_{k}.csv")) for k in SUB) \
        .to_csv(os.path.join(RES, "functional_residues.csv"), index=False)
    print(s.round(4).to_string(index=False))
    for k in SUB:
        t = tables[k]
        c = t[t.pct >= LOW].ligand_state.value_counts().to_dict()
        print(f"   {SUB[k][0]} ligand states of clusters ≥30%: {c}")
    figures(tables, counts)
    print("Tables in", RES, "| figures in", FIG)


if __name__ == "__main__":
    main()
