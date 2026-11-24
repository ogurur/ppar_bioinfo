# =====================================================================
# Which helix is each conserved water next to?
# ---------------------------------------------------------------------
# For every conserved water site (>= MIN_PCT % at 1 Å) of each subtype:
#   - protein residues with any heavy atom within CONTACT Å of the reference
#     water, and the structural element (helix / loop / β-sheet) of each;
#   - the nearest helix and its distance, even if it is farther than CONTACT.
# Helix names use PPARγ1 numbering (1PRG).  PPARα and PPARδ residues are
# translated to PPARγ numbering by sequence alignment, so "H3" means the same
# helix in all three subtypes.
#
# Element boundaries (PPARγ1 numbering): the 1PRG authors' helix records,
# named with literature anchors - H3 276-303, Ω-loop 263-273, H12 466-476,
# Tyr320 on H5, Gln454 on H11.  H6-H9 are named by order: check them
# against a PPARγ structure paper before publishing.
#
# Input  (project folder):
#   data/raw/PPAR<subtype>/<ref>.cif
#   data/results/water_sites_conservation_<subtype>.csv
# Output:
#   data/results/water_site_helices_<subtype>.csv
#   figures/pymol_sites_<subtype>.pml      (run in PyMOL from the project folder)
#   figures/sites_per_element.png
#
# Run:  python label_helices.py       Needs: numpy, pandas, matplotlib, biopython
# =====================================================================
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from Bio import Align

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else BASE_DIR

MIN_PCT = 30.0     # label sites with at least this conservation at 1 Å
CONTACT = 4.0      # Å, residue counts as "next to" the water

SUBTYPES = {
    "gamma": {"label": "PPARγ", "ref": "1prg", "chain": "A", "raw": "data/raw/PPARgamma"},
    "alpha": {"label": "PPARα", "ref": "1k7l", "chain": "A", "raw": "data/raw/PPARalpha"},
    "delta": {"label": "PPARδ", "ref": "1gwx", "chain": "A", "raw": "data/raw/PPARdelta"},
}
GAMMA = "gamma"

# structural elements in PPARγ1 numbering (inclusive ranges), in sequence order
ELEMENTS = [
    ("H1", 207, 225, "helix"), ("loop H1–H2", 226, 229, "loop"),
    ("H2", 230, 237, "helix"), ("loop H2–β1", 238, 246, "loop"),
    ("β-sheet", 247, 249, "sheet"), ("loop β1–H2'", 250, 251, "loop"),
    ("H2'", 252, 259, "helix"), ("Ω-loop", 260, 275, "loop"),
    ("H3", 276, 303, "helix"), ("loop H3–H4/5", 304, 305, "loop"),
    ("H4/5", 306, 333, "helix"), ("loop H5–β", 334, 337, "loop"),
    ("β-sheet", 338, 349, "sheet"), ("H6", 350, 355, "helix"),
    ("loop H6–H7", 356, 358, "loop"), ("H7", 359, 377, "helix"),
    ("loop H7–H8", 378, 380, "loop"), ("H8", 381, 392, "helix"),
    ("loop H8–H9", 393, 402, "loop"), ("H9", 403, 424, "helix"),
    ("loop H9–H10", 425, 430, "loop"), ("H10/11", 431, 459, "helix"),
    ("loop H11–H12", 460, 465, "loop"), ("H12", 466, 477, "helix"),
]
HELIX_COLORS = {  # PyMOL colour names, one per helix
    "H1": "paleyellow", "H2": "grey70", "H2'": "grey60", "H3": "slate", "H4/5": "aquamarine",
    "H6": "wheat", "H7": "lightorange", "H8": "palegreen", "H9": "lightpink",
    "H10/11": "violet", "H12": "salmon",
}

AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
SHORT = {k: k[0] + k[1:].lower() for k in AA}


def element_of(gnum):
    """Element name for a residue number in PPARγ1 numbering."""
    if gnum is None:
        return "unaligned"
    for name, a, b, _ in ELEMENTS:
        if a <= gnum <= b:
            return name
    return "N/C-terminus"


# ---------------------------------------------------------------------
def read_chain(path, chain):
    """Protein heavy atoms of one chain: (resi, resn, atom, xyz) + residue list."""
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
    atoms, residues = [], []
    m0 = None
    for r in rows:
        if len(r) != len(tags):
            continue
        if "pdbx_PDB_model_num" in c:
            m0 = m0 or r[c["pdbx_PDB_model_num"]]
            if r[c["pdbx_PDB_model_num"]] != m0:
                continue
        if r[c["auth_asym_id"]] != chain or r[c.get("auth_comp_id", c["label_comp_id"])] not in AA:
            continue
        if "label_alt_id" in c and r[c["label_alt_id"]] not in (".", "?", "A"):
            continue
        if r[c["type_symbol"]] == "H":
            continue
        resi = int(r[c["auth_seq_id"]])
        resn = r[c.get("auth_comp_id", c["label_comp_id"])]
        name = r[c.get("auth_atom_id", c["label_atom_id"])]
        xyz = [float(r[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]
        atoms.append((resi, resn, name, xyz))
        if name == "CA":
            residues.append((resi, resn))
    return atoms, residues


def to_gamma_numbering(res, res_g):
    """{resi: PPARγ1 residue number} by global sequence alignment."""
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2, -1
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    a = al.align("".join(AA[r[1]] for r in res), "".join(AA[r[1]] for r in res_g))[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(*a.aligned):
        for k in range(a1 - a0):
            m[res[a0 + k][0]] = res_g[b0 + k][0]
    return m


def label_subtype(key, cfg, res_g):
    atoms, res = read_chain(os.path.join(BASE_DIR, cfg["raw"], f"{cfg['ref']}.cif"), cfg["chain"])
    gmap = {r[0]: r[0] for r in res} if key == GAMMA else to_gamma_numbering(res, res_g)
    elem = {r[0]: element_of(gmap.get(r[0])) for r in res}
    xyz = np.array([a[3] for a in atoms])
    a_res = np.array([a[0] for a in atoms])
    a_name = {a[0]: a[1] for a in atoms}

    sites = pd.read_csv(os.path.join(BASE_DIR, "data", "results",
                                     f"water_sites_conservation_{key}.csv"))
    sites = sites[sites["pct_1.0A"] >= MIN_PCT]
    helix_names = [e[0] for e in ELEMENTS if e[3] == "helix"]
    rows = []
    for _, s in sites.iterrows():
        w = np.array([s.ref_x, s.ref_y, s.ref_z])
        d = np.linalg.norm(xyz - w, axis=1)
        # residues within CONTACT, closest first
        near = {}
        for i in np.argsort(d):
            if d[i] > CONTACT:
                break
            near.setdefault(int(a_res[i]), float(d[i]))
        by_elem = {}
        for r, dist in near.items():
            by_elem.setdefault(elem[r], []).append((dist, r))
        elems_sorted = sorted(by_elem.items(), key=lambda kv: min(kv[1])[0])
        contacts = "; ".join(
            f"{e} (" + ", ".join(f"{SHORT[a_name[r]]}{r}" for _, r in sorted(v)) + ")"
            for e, v in elems_sorted)
        # nearest helix, whatever the distance
        best = (None, np.inf)
        for h in set(helix_names):
            idx = np.isin(a_res, [r for r in elem if elem[r] == h])
            if idx.any():
                dh = float(d[idx].min())
                if dh < best[1]:
                    best = (h, dh)
        location = " / ".join(e for e, _ in elems_sorted) or "no residue within %.1f Å" % CONTACT
        rows.append({
            "W": s.W, "ref_water_resi": s.ref_water_resi, "pct_1.0A": s["pct_1.0A"],
            "class": "high" if s["pct_1.0A"] > 70 else "moderate",
            "location": location,
            "nearest_helix": best[0], "dist_nearest_helix_A": round(best[1], 2),
            "residues_within_%gA" % CONTACT: contacts,
            "residues_gamma_numbering": ", ".join(
                f"{SHORT[a_name[r]]}{r}→{gmap.get(r, '?')}" for r in near) if key != GAMMA else "",
        })
    table = pd.DataFrame(rows)
    return table, res, elem, gmap


def pymol_script(key, cfg, table, res, elem):
    """Reference cartoon coloured by helix, helix labels, sites as spheres."""
    path = f"data/raw/{cfg['raw'].split('/')[-1]}/{cfg['ref']}.cif"
    obj, ch = cfg["ref"], cfg["chain"]
    L = [f"# {cfg['label']}: conserved water sites (>= {MIN_PCT:g}% at 1 Å) and nearby helices",
         "# run from the project folder:  pymol figures/pymol_sites_%s.pml" % key,
         "reinitialize", f"load {path}, {obj}", f"remove {obj} and not chain {ch}",
         "bg_color white", "hide everything", f"show cartoon, {obj}",
         f"color grey90, {obj}", "set cartoon_transparency, 0.15",
         "set label_size, 16", "set label_color, black", "set float_labels, 1"]
    for h, col in HELIX_COLORS.items():
        rs = [r for r, _ in res if elem.get(r) == h]
        if not rs:
            continue
        sel = "+".join(str(r) for r in rs)
        tag = h.replace("/", "_").replace("'", "p")
        L += [f"select hx_{tag}, {obj} and resi {sel}", f"color {col}, hx_{tag}",
              f"pseudoatom lab_{tag}, hx_{tag} and name CA", f'label lab_{tag}, "{h}"',
              f"hide everything, lab_{tag}", f"show labels, lab_{tag}"]
    for _, s in table.iterrows():
        col = "red" if s["class"] == "high" else "orange"
        L += [f"select site_{s.W}, {obj} and resn HOH and resi {s.ref_water_resi}",
              f"show spheres, site_{s.W}", f"set sphere_scale, {0.45 if col == 'red' else 0.3}, site_{s.W}",
              f"color {col}, site_{s.W}", f'label site_{s.W}, "{s.W}"']
    L += ["deselect", f"orient {obj}",
          "# when the view is right:  ray 2400, 1800  then  png figures/sites_%s.png, dpi=300" % key]
    return "\n".join(L) + "\n"


def plot_counts(tables, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = [e[0] for e in ELEMENTS if e[3] == "helix"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), dpi=300, sharey=True)
    for ax, (key, t) in zip(axes, tables.items()):
        hi = Counter(t[t["class"] == "high"]["nearest_helix"])
        mo = Counter(t[t["class"] == "moderate"]["nearest_helix"])
        y = np.arange(len(order))
        ax.barh(y, [mo[h] for h in order], color="#9ec5f4", label="moderate (30–70%)")
        ax.barh(y, [hi[h] for h in order], left=[mo[h] for h in order], color="#256abf",
                label="high (>70%)")
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=9)
        ax.invert_yaxis()
        ax.set_title(SUBTYPES[key]["label"], fontsize=11, loc="left")
        ax.set_xlabel("Conserved water sites", fontsize=9)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.xaxis.grid(True, color="#e4e3de")
        ax.set_axisbelow(True)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], fontsize=8, frameon=False, ncol=2,
               loc="upper right", bbox_to_anchor=(0.99, 0.99))
    fig.suptitle("Nearest helix of each conserved water site (≥30% of structures, 1 Å)",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, facecolor="white")


def main():
    res_dir = os.path.join(OUT_DIR, "data", "results")
    fig_dir = os.path.join(OUT_DIR, "figures")
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    g = SUBTYPES[GAMMA]
    _, res_g = read_chain(os.path.join(BASE_DIR, g["raw"], f"{g['ref']}.cif"), g["chain"])
    tables = {}
    for key, cfg in SUBTYPES.items():
        t, res, elem, gmap = label_subtype(key, cfg, res_g)
        t.to_csv(os.path.join(res_dir, f"water_site_helices_{key}.csv"), index=False)
        with open(os.path.join(fig_dir, f"pymol_sites_{key}.pml"), "w", encoding="utf-8") as f:
            f.write(pymol_script(key, cfg, t, res, elem))
        tables[key] = t
        print(f"== {cfg['label']}: {len(t)} sites")
        print(t[["W", "pct_1.0A", "location", "nearest_helix", "dist_nearest_helix_A"]]
              .to_string(index=False))
    plot_counts({k: tables[k] for k in ("alpha", "gamma", "delta")},
                os.path.join(fig_dir, "sites_per_element.png"))
    print("Tables in", res_dir, "| figures and .pml in", fig_dir)


if __name__ == "__main__":
    main()
