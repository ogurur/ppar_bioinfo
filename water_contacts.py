# =====================================================================
# What do the conserved waters touch?  +  Are sites shared by α, γ, δ?
# ---------------------------------------------------------------------
# Part 1 - contacts / bridges (per subtype)
#   For every conserved site (>= MIN_PCT % at 1 Å) and every structure that
#   has a water there, list the polar atoms (N/O) within HB_CUTOFF Å:
#   protein residues, other waters, ligands/ions.  Classify the water:
#     direct bridge        residue > H2O < residue   (H-bonds >= 2 residues)
#     water-mediated       residue > H2O - H2O < residue
#     single residue       H-bonds one residue only
#     waters only / none
#   Residue numbers of every structure are mapped onto the REFERENCE
#   numbering by sequence alignment, so "LYS301" means the same residue in
#   every structure of that subtype.
#
# Part 2 - cross-subtype
#   Superimpose the α (1K7L) and δ (1GWX) references onto γ (1PRG) using
#   Cα of sequence-aligned residues, move their water sites into the γ
#   frame and find sites that coincide (<= SHARED_CUTOFF Å).
#
# Input (project folder):
#   data/raw/PPAR<subtype>/*.cif
#   conserved_align_waters_<subtype>.csv           (1 Å run)
#   data/results/water_sites_conservation_<subtype>.csv  (plot_water_conservation.py)
# Output (data/results/):
#   water_contacts_<subtype>.csv        one row per water occurrence
#   water_site_contacts_<subtype>.csv   one row per conserved site (summary)
#   cross_subtype_sites.csv             sites shared between the 3 references
#   cross_subtype_occupancy.csv         every structure of every subtype moved
#                                       into the γ frame: % with a water at each
#                                       conserved site (reference-independent)
#
# Run:  python water_contacts.py        Needs: numpy, pandas, biopython
# =====================================================================

import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import Align

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR

HB_CUTOFF = 3.5        # Å, polar atom - water oxygen = hydrogen-bond partner
MIN_PCT = 30.0         # analyse sites with at least this conservation at 1 Å
SHARED_CUTOFF = 1.5    # Å, sites of two subtypes closer than this = same site
MIN_SHOW = 20.0        # % - contacts listed in the summary if seen in >= this

SUBTYPES = {
    "gamma": {"label": "PPARγ", "ref": "1prg", "ref_chain": "A", "raw": "data/raw/PPARgamma"},
    "alpha": {"label": "PPARα", "ref": "1k7l", "ref_chain": "A", "raw": "data/raw/PPARalpha"},
    "delta": {"label": "PPARδ", "ref": "1gwx", "ref_chain": "A", "raw": "data/raw/PPARdelta"},
}
CROSS_TARGET = "gamma"   # frame used for the cross-subtype comparison

AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
      "MSE": "M", "SEP": "S", "TPO": "T", "CSO": "C", "CME": "C", "KCX": "K"}
WATER = {"HOH", "WAT", "DOD", "H2O"}
BACKBONE = {"N", "O", "OXT"}


# ---------------------------------------------------------------------
# mmCIF reading (first model, first alt. location)
# ---------------------------------------------------------------------
class Structure:
    def __init__(self, path):
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
        n = len(tags)
        rows = [r for r in rows if len(r) == n]
        if "pdbx_PDB_model_num" in c:
            m0 = rows[0][c["pdbx_PDB_model_num"]]
            rows = [r for r in rows if r[c["pdbx_PDB_model_num"]] == m0]
        alt = c.get("label_alt_id")
        if alt is not None:
            rows = [r for r in rows if r[alt] in (".", "?", "A", "1")]

        def col(name, fallback=None):
            k = c.get(name, c.get(fallback))
            return np.array([r[k] for r in rows])

        self.group = col("group_PDB")
        self.elem = col("type_symbol")
        self.name = col("auth_atom_id", "label_atom_id")
        self.resn = col("auth_comp_id", "label_comp_id")
        self.chain = col("auth_asym_id", "label_asym_id")
        resi = col("auth_seq_id", "label_seq_id")
        if "pdbx_PDB_ins_code" in c:
            ins = col("pdbx_PDB_ins_code")
            resi = np.array([r if i in ("?", ".") else r + i for r, i in zip(resi, ins)])
        self.resi = resi
        self.b = col("B_iso_or_equiv").astype(float)
        self.xyz = np.array([[float(r[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]
                             for r in rows])
        self.is_water = np.isin(self.resn, list(WATER))
        self.is_aa = np.isin(self.resn, list(AA)) & ~self.is_water
        self.polar = np.isin(self.elem, ["N", "O"])

    def residues(self, chain):
        """[(resi, resn, CA xyz)] of amino acids in chain, in file order."""
        sel = np.where(self.is_aa & (self.chain == chain) & (self.name == "CA"))[0]
        return [(self.resi[i], self.resn[i], self.xyz[i]) for i in sel]

    def water_index(self, chain, resi):
        sel = np.where(self.is_water & (self.chain == chain) & (self.resi == str(resi))
                       & (self.elem == "O"))[0]
        return int(sel[0]) if len(sel) else None


_aligner = Align.PairwiseAligner()
_aligner.mode = "global"
_aligner.match_score, _aligner.mismatch_score = 2, -1
_aligner.open_gap_score, _aligner.extend_gap_score = -10, -0.5


def residue_map(res_a, res_b):
    """Sequence-align two residue lists -> {resi_a: index_b}."""
    sa = "".join(AA[r[1]] for r in res_a)
    sb = "".join(AA[r[1]] for r in res_b)
    aln = _aligner.align(sa, sb)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            m[res_a[a0 + k][0]] = b0 + k
    return m


# ---------------------------------------------------------------------
# Part 1 helpers
# ---------------------------------------------------------------------
def nearest_sites(df, ref_sites):
    """Assign each matched water to its nearest reference water (site)."""
    pts = np.array([p for _, ps in ref_sites for p in ps])
    owner = np.array([r for r, ps in ref_sites for _ in ps])
    d = np.linalg.norm(df[["x", "y", "z"]].to_numpy()[:, None] - pts[None], axis=2)
    df = df.copy()
    df["ref_water_resi"] = owner[d.argmin(1)]
    df["dist"] = d.min(1)
    return df


def ref_sites_of(st, chain):
    sites = {}
    for i in np.where(st.is_water & (st.chain == chain) & (st.elem == "O"))[0]:
        sites.setdefault(st.resi[i], []).append(st.xyz[i])
    return list(sites.items())


def polar_partners(st, i_water, label_res):
    """Polar atoms within HB_CUTOFF of atom i_water.
    Returns (protein_atoms, water_idx, other_het) - protein atoms as
    'RES123:ATOM' labels (reference numbering via label_res)."""
    d = np.linalg.norm(st.xyz - st.xyz[i_water], axis=1)
    near = np.where((d <= HB_CUTOFF) & st.polar)[0]
    prot, wat, het = [], [], []
    for j in near:
        if j == i_water:
            continue
        if st.is_water[j]:
            wat.append(j)
        elif st.is_aa[j]:
            prot.append((label_res(st.chain[j], st.resi[j], st.resn[j]), st.name[j]))
        else:
            het.append(st.resn[j])
    return prot, wat, het


def analyse_subtype(key, cfg):
    raw = os.path.join(BASE_DIR, cfg["raw"])
    ref = Structure(os.path.join(raw, f"{cfg['ref']}.cif"))
    ref_res = ref.residues(cfg["ref_chain"])
    sites = ref_sites_of(ref, cfg["ref_chain"])

    table = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                     f"water_sites_conservation_{key}.csv"), dtype={"ref_water_resi": str})
    keep = table[table["pct_1.0A"] >= MIN_PCT]
    w_of = dict(zip(table["ref_water_resi"], table["W"]))
    pct_of = dict(zip(table["ref_water_resi"], table["pct_1.0A"]))

    m = pd.read_csv(os.path.join(BASE_DIR, f"conserved_align_waters_{key}.csv"),
                    dtype={"water_resi": str, "chain": str})
    m = nearest_sites(m, sites)
    m = m[m["ref_water_resi"].isin(keep["ref_water_resi"])]
    # one water per structure per site: the closest one
    m = m.sort_values("dist").drop_duplicates(["pdb_id", "ref_water_resi"])

    # the reference itself also has every one of its waters
    ref_rows = pd.DataFrame({"pdb_id": cfg["ref"], "chain": cfg["ref_chain"],
                             "water_resi": keep["ref_water_resi"].values,
                             "ref_water_resi": keep["ref_water_resi"].values, "dist": 0.0})
    m = pd.concat([ref_rows, m[ref_rows.columns]], ignore_index=True)

    records = []
    for pdb_id, grp in m.groupby("pdb_id"):
        path = os.path.join(raw, f"{pdb_id}.cif")
        if not os.path.isfile(path):
            print(f"  [WARN] {path} missing")
            continue
        st = ref if pdb_id == cfg["ref"] else Structure(path)
        maps = {}

        def label_res(chain, resi, resn, st=st, maps=maps, grp=grp):
            main_chain = grp["chain"].iloc[0]
            if chain != main_chain:
                return f"{resn}{resi}(ch{chain})"          # other molecule in the crystal
            if "m" not in maps:
                maps["m"] = residue_map(st.residues(chain), ref_res)
            k = maps["m"].get(resi)
            return f"{ref_res[k][1]}{ref_res[k][0]}" if k is not None else f"{resn}?{resi}"

        for _, row in grp.iterrows():
            iw = st.water_index(row["chain"], row["water_resi"])
            if iw is None:
                print(f"  [WARN] water {row['water_resi']} chain {row['chain']} not found in {pdb_id}")
                continue
            prot, wat, het = polar_partners(st, iw, label_res)
            residues = sorted({p[0] for p in prot})
            direct = sorted(f"{a}|{b}" for a, b in combinations(residues, 2))
            # residue > H2O - H2O' < residue (residue of H2O' not already bonded)
            via = set()
            for j in wat:
                p2, _, _ = polar_partners(st, j, label_res)
                for r2 in {p[0] for p in p2} - set(residues):
                    for r1 in residues:
                        via.add(f"{r1}|{r2}")
            if len(residues) >= 2:
                pattern = "direct bridge"
            elif via:
                pattern = "water-mediated bridge"
            elif residues:
                pattern = "single residue"
            elif wat:
                pattern = "waters only"
            else:
                pattern = "no polar contact"
            records.append({
                "W": w_of[row["ref_water_resi"]], "ref_water_resi": row["ref_water_resi"],
                "pdb_id": pdb_id, "chain": row["chain"], "water_resi": row["water_resi"],
                "B": st.b[iw], "pattern": pattern,
                "residues": ";".join(residues),
                "atoms": ";".join(sorted(f"{r}:{a}" for r, a in prot)),
                "n_waters_nearby": len(wat),
                "ligand_or_ion": ";".join(sorted(set(het))),
                "direct_bridges": ";".join(direct),
                "water_bridges": ";".join(sorted(via)),
                "has_direct_bridge": len(residues) >= 2,
                "has_water_bridge": bool(via),
            })
    occ = pd.DataFrame(records)

    # ---- per-site summary --------------------------------------------
    def top(counter, n, k=6):
        items = [(x, 100 * c / n) for x, c in counter.most_common() if 100 * c / n >= MIN_SHOW]
        return "; ".join(f"{x} {p:.0f}%" for x, p in items[:k])

    ref_nb = {}
    for r, _ in sites:
        i = ref.water_index(cfg["ref_chain"], r)
        d = np.linalg.norm(ref.xyz[ref.is_water] - ref.xyz[i], axis=1)
        ref_nb[r] = int(((d > 0.01) & (d <= HB_CUTOFF)).sum())

    summ = []
    for r, g in occ.groupby("ref_water_resi"):
        n = len(g)
        cnt = lambda col: Counter(x for v in g[col] for x in v.split(";") if x)
        pat = g["pattern"].value_counts(normalize=True) * 100
        row = table[table["ref_water_resi"] == r].iloc[0]
        summ.append({
            "W": w_of[r], "ref_water_resi": r, "pct_conserved_1A": pct_of[r],
            "pct_0.5A": row.get("pct_0.5A"), "pct_2.0A": row.get("pct_2.0A"),
            "n_structures_analysed": n,
            "residues_contacted": top(cnt("residues"), n),
            "atoms_contacted": top(cnt("atoms"), n),
            "pct_direct_bridge": round(100 * g["has_direct_bridge"].mean(), 1),
            "pct_water_bridge": round(100 * g["has_water_bridge"].mean(), 1),
            "main_pattern": pat.index[0],
            "top_direct_bridge": top(cnt("direct_bridges"), n, 3).replace("|", " > H2O < "),
            "top_water_bridge": top(cnt("water_bridges"), n, 3).replace("|", " > H2O - H2O < "),
            "mean_waters_nearby": round(g["n_waters_nearby"].mean(), 2),
            "waters_nearby_in_ref": ref_nb[r],
            "pct_with_ligand_or_ion": round(100 * (g["ligand_or_ion"] != "").mean(), 1),
            "ligands_ions": top(cnt("ligand_or_ion"), n, 4),
            "mean_B": round(g["B"].mean(), 1),
            "pct_B_over_40": round(100 * (g["B"] > 40).mean(), 1),
        })
    summ = pd.DataFrame(summ).sort_values("pct_conserved_1A", ascending=False)
    return occ, summ, ref, ref_res, sites


# ---------------------------------------------------------------------
# Part 2: cross-subtype
# ---------------------------------------------------------------------
def kabsch(P, Q):
    """R, t minimising |R P + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def superpose(res_mobile, res_target, cycles=5, cutoff=2.0):
    """Cα superposition with PyMOL-like outlier rejection."""
    mp = residue_map(res_mobile, res_target)
    pairs = [(i, mp[r[0]]) for i, r in enumerate(res_mobile) if r[0] in mp]
    P = np.array([res_mobile[i][2] for i, _ in pairs])
    Q = np.array([res_target[j][2] for _, j in pairs])
    keep = np.ones(len(P), bool)
    for _ in range(cycles + 1):
        R, t = kabsch(P[keep], Q[keep])
        d = np.linalg.norm(P @ R.T + t - Q, axis=1)
        new = d <= cutoff
        if new.sum() < 3 or (new == keep).all():
            break
        keep = new
    rmsd = float(np.sqrt((d[keep] ** 2).mean()))
    return R, t, rmsd, int(keep.sum()), len(P), mp


def cross_subtype(results):
    tgt = CROSS_TARGET
    _, s_t, _, res_t, _ = results[tgt]
    frames = {tgt: (np.eye(3), np.zeros(3), None)}
    log = []
    for k in results:
        if k == tgt:
            continue
        R, t, rmsd, n_keep, n_all, mp = superpose(results[k][3], res_t)
        frames[k] = (R, t, mp)
        log.append(f"{SUBTYPES[k]['label']} {SUBTYPES[k]['ref'].upper()} -> "
                   f"{SUBTYPES[tgt]['label']} {SUBTYPES[tgt]['ref'].upper()}: "
                   f"Cα RMSD {rmsd:.2f} Å over {n_keep}/{n_all} residues")

    def to_target_numbering(k, label):
        """'LYS301' in subtype k numbering -> residue in target numbering."""
        if k == tgt or not label or "?" in label or "(ch" in label:
            return label
        mp = frames[k][2]
        resi = label[3:]
        j = mp.get(resi)
        return f"{res_t[j][1]}{res_t[j][0]}" if j is not None else label + "(no eq.)"

    # all analysed sites of every subtype, in the target frame
    pts = []
    for k, (_, summ, ref, _, sites) in results.items():
        R, t, _ = frames[k]
        coords = {r: np.array(ps[0]) for r, ps in sites}
        for _, s in summ.iterrows():
            xyz = R @ coords[s["ref_water_resi"]] + t
            res = [x.rsplit(" ", 1)[0] for x in s["residues_contacted"].split("; ") if x]
            pts.append({"subtype": k, "W": s["W"], "pct": s["pct_conserved_1A"], "xyz": xyz,
                        "contacts": ", ".join(f"{r}" + (f"(={to_target_numbering(k, r)})"
                                                         if k != tgt else "") for r in res[:4])})
    # greedy grouping: sites within SHARED_CUTOFF of each other
    used, rows = set(), []
    order = sorted(range(len(pts)), key=lambda i: -pts[i]["pct"])
    for i in order:
        if i in used:
            continue
        grp = {pts[i]["subtype"]: i}
        for j in order:
            if j in used or j == i or pts[j]["subtype"] in grp:
                continue
            if np.linalg.norm(pts[j]["xyz"] - pts[i]["xyz"]) <= SHARED_CUTOFF:
                grp[pts[j]["subtype"]] = j
        used.update(grp.values())
        row = {"n_subtypes": len(grp), "_anchor": pts[i]["xyz"]}
        for k in results:
            j = grp.get(k)
            lab = SUBTYPES[k]["label"]
            row[f"{lab}_W"] = pts[j]["W"] if j is not None else ""
            row[f"{lab}_pct"] = pts[j]["pct"] if j is not None else np.nan
            row[f"{lab}_contacts"] = pts[j]["contacts"] if j is not None else ""
            row[f"{lab}_dist_A"] = (round(float(np.linalg.norm(pts[j]["xyz"] - pts[i]["xyz"])), 2)
                                    if j is not None else np.nan)
        rows.append(row)
    # For subtypes without a matching conserved site: nearest reference water of
    # that subtype at ANY conservation level (it may exist but be < MIN_PCT, or
    # the reference may simply have no water there).
    all_sites = {}
    for k, (_, _, _, _, sites) in results.items():
        R, t, _ = frames[k]
        tab = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                       f"water_sites_conservation_{k}.csv"),
                          dtype={"ref_water_resi": str})
        info = dict(zip(tab["ref_water_resi"], zip(tab["W"], tab["pct_1.0A"])))
        all_sites[k] = [(info[r][0], info[r][1], R @ np.array(ps[0]) + t) for r, ps in sites]
    for row in rows:
        anchor = row.pop("_anchor")
        for k in results:
            lab = SUBTYPES[k]["label"]
            if row[f"{lab}_W"]:
                row[f"{lab}_nearest_any"] = ""
                continue
            W, pct, d = min(((W, pct, float(np.linalg.norm(xyz - anchor)))
                             for W, pct, xyz in all_sites[k]), key=lambda x: x[2])
            row[f"{lab}_nearest_any"] = f"{W} at {d:.1f} Å ({pct:g}%)"
    df = pd.DataFrame(rows).sort_values(["n_subtypes", f"{SUBTYPES[tgt]['label']}_pct"],
                                        ascending=False)
    return df, log


# ---------------------------------------------------------------------
# Part 3: reference-independent cross-subtype occupancy
#   Every structure of every subtype is superimposed (Cα, sequence-aligned)
#   onto the target reference, and all its LBD waters are moved into that
#   frame.  For each conserved site of any subtype we then count, per
#   subtype, the % of structures with a water within OCC_CUTOFF Å.
# ---------------------------------------------------------------------
OCC_CUTOFF = 1.0


def cross_occupancy(results):
    tgt = CROSS_TARGET
    res_t = results[tgt][3]
    # conserved sites of all subtypes, in the target frame
    sites = []
    for k, (_, summ, ref, ref_res, rsites) in results.items():
        R, t = (np.eye(3), np.zeros(3)) if k == tgt else superpose(ref_res, res_t)[:2]
        coords = {r: np.array(ps[0]) for r, ps in rsites}
        for _, s in summ.iterrows():
            sites.append((f"{SUBTYPES[k]['label']}-{s['W']}", k, s["pct_conserved_1A"],
                          R @ coords[s["ref_water_resi"]] + t))
    S = np.array([x[3] for x in sites])
    hits = {k: np.zeros(len(sites)) for k in results}
    n_struct = {k: 0 for k in results}
    for k, cfg in SUBTYPES.items():
        raw = os.path.join(BASE_DIR, cfg["raw"])
        ch = pd.read_csv(os.path.join(BASE_DIR, f"lbd_chains_PPAR{k}.csv"), dtype=str).fillna("")
        entries = [(r.pdb_id, r.used_chains.split()[0]) for r in ch.itertuples() if r.used_chains]
        entries.append((cfg["ref"], cfg["ref_chain"]))
        for pdb_id, chain in entries:
            path = os.path.join(raw, f"{pdb_id}.cif")
            if not os.path.isfile(path):
                continue
            st = results[k][2] if pdb_id == cfg["ref"] else Structure(path)
            res = st.residues(chain)
            if len(res) < 50:
                continue
            R, t, rmsd, *_ = superpose(res, res_t)
            if rmsd > 2.0:
                continue
            w = st.xyz[st.is_water & (st.chain == chain) & (st.elem == "O")] @ R.T + t
            n_struct[k] += 1
            if len(w):
                d = np.linalg.norm(S[:, None] - w[None], axis=2).min(1)
                hits[k] += d <= OCC_CUTOFF
    rows = []
    for i, (name, k, pct, _) in enumerate(sites):
        row = {"site": name, "defined_in": SUBTYPES[k]["label"], "pct_in_own_analysis": pct}
        for k2 in results:
            row[f"pct_in_{SUBTYPES[k2]['label']}"] = round(100 * hits[k2][i] / n_struct[k2], 1)
        rows.append(row)
    df = pd.DataFrame(rows)
    pc = [c for c in df.columns if c.startswith("pct_in_PPAR")]
    df["n_subtypes_over_50pct"] = (df[pc] >= 50).sum(axis=1)
    df["min_pct_all_subtypes"] = df[pc].min(axis=1)
    df["x"], df["y"], df["z"] = S[:, 0].round(2), S[:, 1].round(2), S[:, 2].round(2)
    # the same spot found from two subtypes' references -> mark the duplicate
    df["mean_pct"] = df[pc].mean(axis=1)
    df = df.sort_values(["mean_pct"], ascending=False).reset_index(drop=True)
    same = [""] * len(df)
    xyz = df[["x", "y", "z"]].to_numpy()
    for i in range(len(df)):
        for j in range(i):
            if not same[j] and np.linalg.norm(xyz[i] - xyz[j]) <= SHARED_CUTOFF:
                same[i] = df.loc[j, "site"]
                break
    df["same_site_as"] = same
    df = df.drop(columns="mean_pct")
    return df, n_struct


# ---------------------------------------------------------------------
def main():
    res_dir = os.path.join(OUT_DIR, "data", "results")
    os.makedirs(res_dir, exist_ok=True)
    results = {}
    for k, cfg in SUBTYPES.items():
        print(f"== {cfg['label']}")
        occ, summ, ref, ref_res, sites = analyse_subtype(k, cfg)
        occ.to_csv(os.path.join(res_dir, f"water_contacts_{k}.csv"), index=False)
        summ.to_csv(os.path.join(res_dir, f"water_site_contacts_{k}.csv"), index=False)
        print(f"   {len(summ)} sites >= {MIN_PCT:g}%, {len(occ)} water occurrences analysed")
        results[k] = (occ, summ, ref, ref_res, sites)

    cross, log = cross_subtype(results)
    cross.to_csv(os.path.join(res_dir, "cross_subtype_sites.csv"), index=False)
    print("\n".join(log))
    print(f"   shared by 3 subtypes: {(cross.n_subtypes == 3).sum()}, "
          f"by 2: {(cross.n_subtypes == 2).sum()}")
    occ3, n_struct = cross_occupancy(results)
    occ3.to_csv(os.path.join(res_dir, "cross_subtype_occupancy.csv"), index=False)
    print("   structures superimposed on", SUBTYPES[CROSS_TARGET]["ref"].upper(), n_struct)
    uniq = occ3[occ3.same_site_as == ""]
    print(f"   unique sites: {len(uniq)}; with >= 50% in all three subtypes: "
          f"{(uniq.n_subtypes_over_50pct == 3).sum()}")
    print(f"Tables in {res_dir}")


if __name__ == "__main__":
    main()
