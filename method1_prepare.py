# =====================================================================
# Method 1, step 1: superimpose every LBD chain and extract ALL its waters
# ---------------------------------------------------------------------
# Superposition (project description, section 4.2):
#   Cα atoms of the core LBD helices only - H3, H4/5, H7, H10/11 - no
#   outlier rejection, so the RMSD is the honest core RMSD.  Helices are
#   defined in PPARγ1 numbering (1PRG) and carried to PPARα / PPARδ by
#   sequence alignment, so the same helices are used for every subtype.
#   For the sensitivity analysis the full-LBD Cα superposition is also stored.
# Waters (section 4.3): every water oxygen of the LBD chain within SHELL Å of
#   that chain's protein atoms; original coordinates + B-factor are saved,
#   superpositions are saved as rotation/translation per structure, so the
#   clustering step can place waters in any reference frame.
# References: the main one per subtype (1K7L, 1PRG, 1GWX) plus N_ALT
#   alternatives (highest resolution, >= MIN_WATERS waters, complete chain)
#   for the reference-choice sensitivity test.
#
# Input (project folder): data/raw/PPAR<subtype>/*.cif, lbd_chains_PPAR<subtype>.csv,
#                         data/metadata.csv
# Output (data/processed/method1/):
#   waters_<subtype>.csv      one row per water, original coordinates
#   transforms_<subtype>.csv  one row per structure x reference x scope (R, t, RMSD)
#   references.csv            main and alternative references
#
# Run:  python method1_prepare.py      Needs: numpy, pandas, scipy, biopython
# =====================================================================
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from Bio import Align

BASE_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE_DIR, "data", "processed", "method1")

SHELL = 5.0          # Å, keep waters within this distance of the chain's protein atoms
N_ALT = 2            # alternative references per subtype
MIN_WATERS = 50      # an alternative reference must have at least this many LBD waters
MIN_RES_LEN = 250    # ... and at least this many residues with Cα

SUBTYPES = {
    "alpha": {"label": "PPARα", "ref": "1k7l", "chain": "A", "raw": "data/raw/PPARalpha", "meta": "PPARalpha"},
    "gamma": {"label": "PPARγ", "ref": "1prg", "chain": "A", "raw": "data/raw/PPARgamma", "meta": "PPARgamma"},
    "delta": {"label": "PPARδ", "ref": "1gwx", "chain": "A", "raw": "data/raw/PPARdelta", "meta": "PPARdelta"},
}
GAMMA = "gamma"
# core helices in PPARγ1 numbering (1PRG helix records; H3 276-303 per literature)
CORE = [(276, 303), (306, 333), (359, 377), (431, 459)]     # H3, H4/5, H7, H10/11

AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
WATER = {"HOH", "WAT", "DOD"}


def read_chain(path, chain):
    """Cα list, protein heavy-atom coords and water oxygens of one chain."""
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
    rn = c.get("auth_comp_id", c["label_comp_id"])
    an = c.get("auth_atom_id", c["label_atom_id"])
    ca, prot, wat, m0 = [], [], [], None
    seen_w = set()
    for r in rows:
        if len(r) != len(tags) or r[c["auth_asym_id"]] != chain:
            continue
        if "pdbx_PDB_model_num" in c:
            m0 = m0 or r[c["pdbx_PDB_model_num"]]
            if r[c["pdbx_PDB_model_num"]] != m0:
                continue
        if "label_alt_id" in c and r[c["label_alt_id"]] not in (".", "?", "A", "1"):
            continue
        if r[c["type_symbol"]] in ("H", "D"):
            continue
        xyz = [float(r[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]
        b = float(r[c["B_iso_or_equiv"]])
        resn, resi = r[rn], r[c["auth_seq_id"]]
        if resn in AA:
            prot.append(xyz)
            if r[an] == "CA":
                ca.append((resi, resn, xyz, b))
        elif resn in WATER and r[c["type_symbol"]] == "O" and resi not in seen_w:
            seen_w.add(resi)
            occ = float(r[c["occupancy"]]) if "occupancy" in c else 1.0
            wat.append((resi, xyz, b, occ))
    return ca, np.array(prot), wat


_al = Align.PairwiseAligner()
_al.mode = "global"
_al.match_score, _al.mismatch_score = 2, -1
_al.open_gap_score, _al.extend_gap_score = -10, -0.5


def seq_map(a, b):
    """index map a -> b of two Cα lists by global sequence alignment."""
    aln = _al.align("".join(AA[x[1]] for x in a), "".join(AA[x[1]] for x in b))[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            m[a0 + k] = b0 + k
    return m


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def superpose(ca, ref_ca, ref_core_idx):
    """(core R, t, rmsd, n), (full R, t, rmsd, n) onto the reference."""
    m = seq_map(ca, ref_ca)
    out = []
    for scope in ("core", "full"):
        pairs = [(i, j) for i, j in m.items() if scope == "full" or j in ref_core_idx]
        if len(pairs) < 20:
            out.append((None, None, np.nan, len(pairs)))
            continue
        P = np.array([ca[i][2] for i, _ in pairs])
        Q = np.array([ref_ca[j][2] for _, j in pairs])
        R, t, rmsd = kabsch(P, Q)
        out.append((R, t, rmsd, len(pairs)))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = pd.read_csv(os.path.join(BASE_DIR, "data", "metadata.csv"))
    meta["pdb_id"] = meta["pdb_id"].str.lower()
    meta = meta.set_index("pdb_id")

    # PPARγ reference residues -> core definition carried to other subtypes
    g = SUBTYPES[GAMMA]
    gca, _, _ = read_chain(os.path.join(BASE_DIR, g["raw"], f"{g['ref']}.cif"), g["chain"])

    ref_rows = []
    for key, cfg in SUBTYPES.items():
        raw = os.path.join(BASE_DIR, cfg["raw"])
        ch = pd.read_csv(os.path.join(BASE_DIR, f"lbd_chains_PPAR{key}.csv"), dtype=str).fillna("")
        entries = [(r.pdb_id.lower(), r.used_chains.split()[0]) for r in ch.itertuples() if r.used_chains]
        entries = [(cfg["ref"], cfg["chain"])] + [e for e in entries if e[0] != cfg["ref"]]

        # parse every structure once
        parsed = {}
        for pid, chain in entries:
            path = os.path.join(raw, f"{pid}.cif")
            if os.path.isfile(path):
                parsed[pid] = (chain,) + read_chain(path, chain)
        print(f"== {cfg['label']}: parsed {len(parsed)} structures")

        # alternative references: best resolution, enough waters, complete chain
        cand = []
        for pid, (chain, ca, prot, wat) in parsed.items():
            if pid == cfg["ref"] or len(ca) < MIN_RES_LEN or len(wat) < MIN_WATERS:
                continue
            res = meta.loc[pid, "resolution_A"] if pid in meta.index else np.nan
            cand.append((res, pid))
        alts = [p for _, p in sorted(cand)[:N_ALT]]
        refs = [cfg["ref"]] + alts
        for i, r in enumerate(refs):
            ref_rows.append({"subtype": key, "reference": r, "role": "main" if i == 0 else "alternative",
                             "chain": parsed[r][0],
                             "resolution_A": meta.loc[r, "resolution_A"] if r in meta.index else np.nan,
                             "n_waters": len(parsed[r][3]), "n_residues": len(parsed[r][1])})

        # core residues of each reference, via PPARγ numbering
        core_idx = {}
        for r in refs:
            rca = parsed[r][1]
            mg = seq_map(rca, gca)
            core_idx[r] = {i for i, j in mg.items()
                           if any(a <= int(gca[j][0]) <= b for a, b in CORE)}

        # waters (original frame) and transforms onto every reference
        wrows, trows = [], []
        for pid, (chain, ca, prot, wat) in parsed.items():
            if not len(ca):
                continue
            ca_mean = float(np.mean([x[3] for x in ca]))
            if wat and len(prot):
                tree = cKDTree(prot)
                wxyz = np.array([w[1] for w in wat])
                dmin, _ = tree.query(wxyz, k=1)
                for w, dd in zip(wat, dmin):
                    if dd <= SHELL:
                        wrows.append({"pdb_id": pid, "chain": chain, "water_resi": w[0],
                                      "x": w[1][0], "y": w[1][1], "z": w[1][2], "B": w[2],
                                      "occupancy": w[3], "B_norm": round(w[2] / ca_mean, 3),
                                      "dist_protein": round(float(dd), 2)})
            for r in refs:
                (Rc, tc, rc, nc), (Rf, tf, rf, nf) = superpose(ca, parsed[r][1], core_idx[r])
                for scope, R, t, rmsd, n in (("core", Rc, tc, rc, nc), ("full", Rf, tf, rf, nf)):
                    row = {"pdb_id": pid, "reference": r, "scope": scope, "rmsd": round(rmsd, 3),
                           "n_pairs": n}
                    if R is not None:
                        row.update({f"R{i}{j}": R[i, j] for i in range(3) for j in range(3)})
                        row.update({"t0": t[0], "t1": t[1], "t2": t[2]})
                    trows.append(row)
        W = pd.DataFrame(wrows)
        T = pd.DataFrame(trows)
        W["resolution_A"] = W["pdb_id"].map(lambda p: meta.loc[p, "resolution_A"] if p in meta.index else np.nan)
        W.to_csv(os.path.join(OUT, f"waters_{key}.csv"), index=False)
        T.to_csv(os.path.join(OUT, f"transforms_{key}.csv"), index=False)
        main_core = T[(T.reference == cfg["ref"]) & (T.scope == "core") & (T.pdb_id != cfg["ref"])]
        main_full = T[(T.reference == cfg["ref"]) & (T.scope == "full") & (T.pdb_id != cfg["ref"])]
        print(f"   references: {refs}")
        print(f"   waters kept (≤{SHELL:g} Å of protein): {len(W)} in {W.pdb_id.nunique()} structures")
        print(f"   core Cα RMSD to {cfg['ref'].upper()}: median {main_core.rmsd.median():.2f} Å, "
              f"max {main_core.rmsd.max():.2f}, >2 Å: {(main_core.rmsd > 2).sum()} "
              f"(median {main_core.n_pairs.median():.0f} Cα pairs)")
        print(f"   full Cα RMSD: median {main_full.rmsd.median():.2f} Å, max {main_full.rmsd.max():.2f}")
    pd.DataFrame(ref_rows).to_csv(os.path.join(OUT, "references.csv"), index=False)
    print("Saved to", OUT)


if __name__ == "__main__":
    main()
