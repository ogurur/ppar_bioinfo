# =====================================================================
# PPAR water alignment - LBD chains only
# ---------------------------------------------------------------------
# One script for all three subtypes (replaces alignment_pparalpha.py,
# alignment_ppargamma.py, alignment_ppardelta.py).
#
# What changed vs. the old scripts:
#   * Only PPAR ligand-binding-domain (LBD) chains are aligned.  RXR,
#     coactivator / corepressor peptides etc. are skipped.
#   * LBD chains are identified from the mmCIF header (UniProt accession
#     of each entity).  If the header can't be read, falls back to a
#     length check (>= MIN_LBD_RESIDUES CA atoms).
#   * The reference is also restricted to its LBD chain (and its waters).
#   * Entries with several PPAR chains (e.g. A + B) are handled by
#     CHAIN_MODE - default "best" keeps one chain per PDB entry, so a
#     structure is never counted twice.
#   * Writes lbd_chains_<subtype>.csv: which chains were found / used
#     for every entry (useful for metadata.csv).
#   * Superposition uses ONLY cmd.align (sequence alignment first, then
#     fit of the paired residues).  cmd.super is not used anywhere.
#
# How to run (in PyMOL):   run alignment_ppar_lbd.py
# No extra packages needed - only PyMOL's own Python.
# =====================================================================

import os
import csv
import gzip
from pymol import cmd

# ---------------------------------------------------------------------
# Config - edit these for your run
# ---------------------------------------------------------------------
# Which subtypes to run in this go, any of: "PPARalpha", "PPARgamma", "PPARdelta"
SUBTYPES_TO_RUN = ["PPARalpha", "PPARgamma", "PPARdelta"]

# Folder that contains PPARalpha/, PPARgamma/, PPARdelta/ with <id>.cif files.
# Set to None (or "") to download with cmd.fetch instead.
LOCAL_CIF_ROOT = "C:/Users/Ogurur/ppar_bioinfo/data/raw"

CUTOFF = 0.5              # A - waters within this distance of a reference water
CHAIN_MODE = "best"       # "best"  = one LBD chain per entry, lowest align RMSD to reference
                          # "first" = one LBD chain per entry, first in the file
                          # "all"   = every LBD chain (entries with A+B counted twice!)
MIN_LBD_RESIDUES = 150    # fallback length check: LBD ~250-290 res, peptides ~10-30

PPAR_UNIPROT = {
    "PPARalpha": "Q07869",
    "PPARgamma": "P37231",
    "PPARdelta": "Q03181",
}

CONFIG = {
    "PPARalpha": {
        "reference_pdb": "1k7l",
        "reference_chain": None,   # None = first PPAR LBD chain found in the reference
        "output_csv": "conserved_align_waters_alpha_05.csv",
        "log_txt": "conserved_align_waters_log_alpha_05.txt",
        "pdb_ids": [
            "1i7g", "1kkq", "2npa", "2p54", "2rew", "2znn", "3et1",
            "3fei", "3g8i", "3kdt", "3kdu", "3sp6", "3vi8", "4bcr", "4ci4",
            "5azt", "5hyk", "6kax", "6kay", "6kaz", "6kb0", "6kb1", "6kb2",
            "6kb3", "6kb4", "6kb5", "6kb6", "6kb7", "6kb8", "6kb9", "6kba",
            "6kxx", "6kxy", "6kyp", "6l36", "6l37", "6l38", "6l96", "6lx4",
            "6lx5", "6lx6", "6lx7", "6lx8", "6lx9", "6lxa", "6lxb", "6lxc",
            "7bpy", "7bpz", "7bq0", "7bq1", "7bq2", "7bq3", "7bq4", "7c6q",
            "7e5f", "7e5g", "7e5h", "7e5i", "8huk", "8hun", "8huq", "8rce",
            "8yt6", "8yt9", "8ytl", "8yu8", "8ywu", "8ywv", "9iwj", "9iwl",
            "9iwm", "9iwn", "9iwo", "9sfs", "9t5s", "9vzs", "9vzt",
        ],
    },
    "PPARgamma": {
        "reference_pdb": "1prg",   # or "2prg"
        "reference_chain": "A",
        "output_csv": "conserved_align_waters_gamma_05.csv",
        "log_txt": "conserved_align_waters_log_gamma_05.txt",
        "pdb_ids": [
            "1fm6", "1fm9", "1i7i", "1k74", "1knu", "1nyx", "1prg", "1rdt",
            "1wm0", "1zeo", "1zgy", "2ath", "2f4b", "2fvj", "2g0g", "2g0h",
            "2gtk", "2hfp", "2hwq", "2hwr", "2i4j", "2i4p", "2i4z", "2om9",
            "2p4y", "2pob", "2prg", "2q59", "2q5p", "2q5s", "2q61", "2q6r",
            "2q6s", "2q8s", "2vsr", "2vst", "2vv0", "2vv1", "2vv2", "2vv3",
            "2vv4", "2xkw", "2yfe", "2zk0", "2zk1", "2zk2", "2zk3", "2zk4",
            "2zk5", "2zk6", "2zno", "2zvt", "3ads", "3adt", "3adu", "3adv",
            "3adw", "3adx", "3an3", "3an4", "3b0q", "3b0r", "3b1m", "3b3k",
            "3bc5", "3cdp", "3cds", "3cs8", "3cwd", "3d6d", "3dzu", "3dzy",
            "3e00", "3et0", "3et3", "3fej", "3fur", "3g9e", "3gbk", "3h0a",
            "3ho0", "3hod", "3ia6", "3k8s", "3kmg", "3lmp", "3noa", "3osi",
            "3osw", "3pba", "3po9", "3prg", "3qt0", "3r5n", "3r8a", "3r8i",
            "3s9s", "3sz1", "3t03", "3ty0", "3u9q", "3v9t", "3v9v", "3v9y",
            "3vjh", "3vji", "3vn2", "3vso", "3vsp", "3wj4", "3wj5", "3wmh",
            "3x1h", "3x1i", "4a4v", "4a4w", "4ci5", "4e4k", "4e4q", "4em9",
            "4ema", "4f9m", "4fgy", "4hee", "4jaz", "4jl4", "4l96", "4l98",
            "4o8f", "4oj4", "4prg", "4pvu", "4pwl", "4r06", "4r2u", "4r6s",
            "4xld", "4xta", "4xuh", "4xum", "4y29", "4yt1", "5azv", "5dsh",
            "5dv3", "5dv6", "5dv8", "5dvc", "5dwl", "5f9b", "5gtn", "5gto",
            "5gtp", "5hzc", "5ji0", "5lsg", "5tto", "5two", "5u5l", "5ugm",
            "5wqx", "5wr0", "5wr1", "5y2o", "5y2t", "5ycn", "5ycp", "5z5s",
            "5z6s", "6ad9", "6an1", "6aug", "6avi", "6c1i", "6c5q", "6c5t",
            "6d3e", "6d8x", "6d94", "6dbh", "6dcu", "6dgl", "6dgo", "6dgp",
            "6dgq", "6dgr", "6dh9", "6dha", "6e5a", "6enq", "6f2l", "6fzf",
            "6fzg", "6fzj", "6fzp", "6fzy", "6icj", "6ijr", "6ijs", "6ilq",
            "6izm", "6izn", "6jey", "6jf0", "6jq7", "6k0t", "6ktm", "6ktn",
            "6l89", "6l8b", "6mcz", "6md0", "6md1", "6md2", "6md4", "6ms7",
            "6o67", "6o68", "6oni", "6onj", "6pdz", "6qj5", "6t1s", "6t1v",
            "6t6b", "6t9c", "6tdc", "6tsg", "6vzl", "6vzm", "6vzn", "6vzo",
            "6y3u", "6zly", "7a7h", "7ahj", "7awc", "7awd", "7cxe", "7cxf",
            "7cxg", "7cxh", "7cxi", "7cxj", "7cxk", "7cxl", "7e0a", "7e2o",
            "7efq", "7jqg", "7lot", "7p4e", "7qb1", "7rle", "7sqa", "7sqb",
            "7wgo", "7wgp", "7wgq", "7wox", "8adf", "8aqm", "8aqn", "8aty",
            "8atz", "8b8w", "8b8x", "8b8y", "8b8z", "8b90", "8b91", "8b92",
            "8b93", "8b94", "8b95", "8bf1", "8bf2", "8bff", "8c0c", "8cph",
            "8cpi", "8cpj", "8dk4", "8dkn", "8dkv", "8dsy", "8dsz", "8fhe",
            "8fhf", "8fhg", "8fkc", "8fkd", "8fke", "8fkf", "8fkg", "8hhp",
            "8hhq", "8hum", "8hup", "8pbo", "8rej", "8sc9", "8u57", "8wfe",
            "8zfn", "8zfo", "8zfp", "8zfq", "8zfr", "8zfs", "8zft", "9ck0",
            "9cwn", "9f7w", "9f7x", "9gfu", "9gvp", "9gvs", "9gvt", "9gw3",
            "9gw4", "9gw7", "9gw8", "9gwb", "9gwc", "9gwe", "9gwf", "9gwg",
            "9gwh", "9gwi", "9gwk", "9gwl", "9gwm", "9gwn", "9gwp", "9gwq",
            "9gwr", "9gws", "9gwx", "9gwy", "9gx0", "9gx1", "9gx2", "9gx3",
            "9gx4", "9gx5", "9gx6", "9gx7", "9gxd", "9gz7", "9gz8", "9gz9",
            "9gza", "9gzb", "9gzc", "9gze", "9gzf", "9gzi", "9hx2", "9iwk",
            "9l3t", "9mah", "9o9n", "9olc", "9qfz", "9qg0", "9qg2", "9qg3",
            "9qg4", "9qg5", "9qg6", "9qgj", "9r6i", "9r6k", "9v8d", "9v8e",
            "9v8f", "9v8g", "9v8h",
        ],
    },
    "PPARdelta": {
        "reference_pdb": "1gwx",
        "reference_chain": None,
        "output_csv": "conserved_align_waters_delta_05.csv",
        "log_txt": "conserved_align_waters_log_delta_05.txt",
        "pdb_ids": [
            "1y0s", "2awh", "2b50", "2baw", "2gwx", "2j14", "2q5g",
            "2xyj", "2xyw", "2xyx", "2znp", "2znq", "3d5f", "3dy6", "3et2",
            "3gwx", "3gz9", "3oz0", "3peq", "3sp9", "3tkm", "5u3q", "5u3r",
            "5u3s", "5u3t", "5u3u", "5u3v", "5u3w", "5u3x", "5u3y", "5u3z",
            "5u40", "5u41", "5u42", "5u43", "5u44", "5u45", "5u46", "5xmx",
            "5y7x", "5zxi", "6a6p", "7f80", "7vwe", "7vwf", "7vwg", "7vwh",
            "7w0g", "7wgl", "7wgn", "8hf8", "8hul", "8huo", "8znn",
        ],
    },
}


# =====================================================================
# Logging
# =====================================================================
_log_file = None


def log(message=""):
    """Print to the console AND write to the log file."""
    print(message)
    if _log_file is not None:
        _log_file.write(str(message) + "\n")
        _log_file.flush()


# =====================================================================
# Minimal mmCIF header reader (no Biopython needed inside PyMOL)
# =====================================================================
def _cif_tokens(text):
    """Yield (value, is_quoted) tokens from mmCIF text.
    Handles 'quoted', "quoted" and ;multi-line; values."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            yield ("\n".join(buf).strip(), True)
            i += 1
            continue
        j, n = 0, len(line)
        while j < n:
            c = line[j]
            if c.isspace():
                j += 1
                continue
            if c == "#":
                break
            if c in "'\"":
                k = j + 1
                while k < n and not (line[k] == c and (k + 1 == n or line[k + 1].isspace())):
                    k += 1
                yield (line[j + 1:k], True)
                j = k + 1
            else:
                k = j
                while k < n and not line[k].isspace():
                    k += 1
                yield (line[j:k], False)
                j = k
        i += 1


def read_cif_categories(path, wanted=("_struct_ref", "_entity_poly", "_entity")):
    """Return {tag: [values]} for the wanted categories.
    Only reads the header part (stops before _atom_site) so it is fast."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        text = f.read()
    cut = text.find("_atom_site.")
    if cut > 0:
        text = text[:cut]

    toks = list(_cif_tokens(text))
    data = {}
    i = 0
    while i < len(toks):
        val, quoted = toks[i]
        if not quoted and val.lower() == "loop_":
            i += 1
            tags = []
            while i < len(toks) and not toks[i][1] and toks[i][0].startswith("_"):
                tags.append(toks[i][0])
                i += 1
            vals = []
            while i < len(toks):
                v, q = toks[i]
                if not q and (v.startswith("_") or v.lower() == "loop_"
                              or v.lower().startswith("data_")):
                    break
                vals.append(v)
                i += 1
            if tags and any(t.split(".")[0] in wanted for t in tags):
                for k, tag in enumerate(tags):
                    data[tag] = vals[k::len(tags)]
        elif not quoted and val.startswith("_"):
            if i + 1 < len(toks) and val.split(".")[0] in wanted:
                data[val] = [toks[i + 1][0]]
            i += 2
        else:
            i += 1
    return data


def lbd_chains_from_header(cif_path, subtype):
    """Author chain IDs whose entity maps to the PPAR UniProt accession.
    Returns (ppar_chains, description_by_chain)."""
    d = read_cif_categories(cif_path)
    target = PPAR_UNIPROT[subtype]

    ppar_entities = {
        e for e, a in zip(d.get("_struct_ref.entity_id", []),
                          d.get("_struct_ref.pdbx_db_accession", []))
        if a.upper() == target
    }
    names = dict(zip(d.get("_entity.id", []),
                     d.get("_entity.pdbx_description", [])))

    ppar_chains, desc = [], {}
    for e, ch in zip(d.get("_entity_poly.entity_id", []),
                     d.get("_entity_poly.pdbx_strand_id", [])):
        for c in ch.split(","):
            c = c.strip()
            desc[c] = names.get(e, "?")
            if e in ppar_entities:
                ppar_chains.append(c)
    return ppar_chains, desc


# =====================================================================
# Loading and chain selection
# =====================================================================
def find_cif_file(pdb_id, subtype):
    """Locate the mmCIF on disk (local folder or wherever cmd.fetch saved it)."""
    folders = []
    if LOCAL_CIF_ROOT:
        folders.append(os.path.join(LOCAL_CIF_ROOT, subtype))
    try:
        folders.append(cmd.get("fetch_path"))
    except Exception:
        pass
    folders.append(os.getcwd())
    for folder in folders:
        for name in (f"{pdb_id}.cif", f"{pdb_id.upper()}.cif",
                     f"{pdb_id}.cif.gz", f"{pdb_id.upper()}.cif.gz"):
            path = os.path.join(folder, name)
            if os.path.isfile(path):
                return path
    return None


def load_structure(pdb_id, subtype):
    if LOCAL_CIF_ROOT:
        path = find_cif_file(pdb_id, subtype)
        if path is None:
            raise FileNotFoundError(
                f"{pdb_id}.cif not found in {os.path.join(LOCAL_CIF_ROOT, subtype)}")
        cmd.load(path, pdb_id)
    else:
        cmd.fetch(pdb_id, type="cif")


def n_ca(obj, chain):
    return cmd.count_atoms(f"{obj} and chain {chain} and polymer.protein and name CA")


def identify_lbd_chains(pdb_id, subtype):
    """Return (lbd_chains, method, note).  method = accession / length / none."""
    protein_chains = [c for c in cmd.get_chains(pdb_id)
                      if cmd.count_atoms(f"{pdb_id} and chain {c} and polymer.protein") > 0]

    cif_path = find_cif_file(pdb_id, subtype)
    if cif_path:
        try:
            header_chains, desc = lbd_chains_from_header(cif_path, subtype)
            chains = [c for c in header_chains if c in protein_chains]
            if chains:
                short = [c for c in chains if n_ca(pdb_id, c) < MIN_LBD_RESIDUES]
                note = ""
                if short:
                    note = f"PPAR chain(s) {short} shorter than {MIN_LBD_RESIDUES} res - check"
                skipped = [f"{c}={desc.get(c, '?')}" for c in protein_chains
                           if c not in chains]
                if skipped:
                    note = (note + "; " if note else "") + "skipped: " + ", ".join(skipped)
                return chains, "accession", note
        except Exception as exc:
            log(f"  [WARN] could not read header of {cif_path}: {exc}")

    # Fallback: length check (cannot tell RXR from PPAR!)
    chains = [c for c in protein_chains if n_ca(pdb_id, c) >= MIN_LBD_RESIDUES]
    if chains:
        return chains, "length", "header not usable - length fallback, verify manually"
    return [], "none", "no LBD-sized chain found"


def chain_sel(obj, chain):
    return f"{obj} and chain {chain}" if chain else obj


def align_rmsd(mobile, target, transform=0):
    """Sequence-based superposition with cmd.align.
    transform=0 -> only measure, transform=1 -> also move `mobile`.
    Returns (rmsd_after_refinement, n_atoms_after_refinement,
             rmsd_before_refinement, n_atoms_before_refinement)."""
    r = cmd.align(f"{mobile} and polymer.protein", f"{target} and polymer.protein",
                  transform=transform)
    return r[0], r[1], r[3], r[4]


# =====================================================================
# Main per-chain work
# =====================================================================
def align_and_collect(pdb_id, chain, reference, ref_waters_sel, rows):
    """Copy one LBD chain (with its waters), superimpose it on the reference
    LBD and record waters within CUTOFF of a reference water."""
    work_obj = f"{pdb_id}_{chain or 'x'}"
    cmd.create(work_obj, chain_sel(pdb_id, chain))

    if cmd.count_atoms(f"{work_obj} and polymer.protein") == 0:
        log(f"  [skip] {work_obj}: no protein atoms")
        cmd.delete(work_obj)
        return None

    try:
        # one call: aligns sequences, fits paired residues and moves the
        # whole object (protein + its waters) onto the reference
        rmsd, n_aligned, rmsd_raw, n_raw = align_rmsd(work_obj, reference, transform=1)
    except Exception as exc:
        log(f"  [ERROR] {work_obj} alignment failed: {exc}")
        cmd.delete(work_obj)
        return None

    sel_name = f"matched_{work_obj}"
    cmd.select(sel_name, f"({work_obj} and resn HOH) within {CUTOFF} of {ref_waters_sel}")
    model = cmd.get_model(sel_name)
    for atom in model.atom:
        rows.append({
            "pdb_id": pdb_id,
            "chain": chain or "single",
            "water_resi": atom.resi,
            "x": round(atom.coord[0], 2),
            "y": round(atom.coord[1], 2),
            "z": round(atom.coord[2], 2),
            "align_rmsd": round(rmsd, 3),
            "n_aligned_atoms": n_aligned,
            "align_rmsd_before_refine": round(rmsd_raw, 3),
            "n_atoms_before_refine": n_raw,
        })

    n_waters = cmd.count_atoms(f"{work_obj} and resn HOH")
    log(f"  chain {chain}: {len(model.atom)} matched / {n_waters} waters "
        f"(align RMSD={rmsd:.3f} over {n_aligned} atoms; "
        f"before refinement {rmsd_raw:.3f} over {n_raw})")
    if rmsd > 2.0:
        log(f"  [OUTLIER] {pdb_id} chain {chain}: RMSD > 2.0 A - inspect manually")

    cmd.delete(work_obj)
    cmd.delete(sel_name)
    return rmsd


def pick_chains(pdb_id, lbd_chains, reference):
    """Apply CHAIN_MODE to entries that have more than one LBD chain."""
    if CHAIN_MODE == "all" or len(lbd_chains) <= 1:
        return lbd_chains
    if CHAIN_MODE == "first":
        return lbd_chains[:1]
    # "best": lowest align RMSD to reference
    scores = []
    for c in lbd_chains:
        tmp = f"_tmp_{pdb_id}_{c}"
        cmd.create(tmp, chain_sel(pdb_id, c))
        try:
            scores.append((align_rmsd(tmp, reference)[0], c))
        except Exception:
            pass
        cmd.delete(tmp)
    if not scores:
        return lbd_chains[:1]
    scores.sort()
    return [scores[0][1]]


def run_subtype(subtype):
    global _log_file
    cfg = CONFIG[subtype]
    ref_pdb = cfg["reference_pdb"].lower()
    log_path = os.path.abspath(cfg["log_txt"])
    csv_path = os.path.abspath(cfg["output_csv"])
    chains_csv = os.path.abspath(f"lbd_chains_{subtype}.csv")
    _log_file = open(log_path, "w")

    cmd.reinitialize()
    log(f"===== {subtype} =====")
    log(f"Working directory: {os.getcwd()}")
    log(f"CHAIN_MODE={CHAIN_MODE}  alignment=cmd.align (sequence)  CUTOFF={CUTOFF} A")

    # ---- reference: LBD chain only ----------------------------------
    load_structure(ref_pdb, subtype)
    ref_lbd, method, note = identify_lbd_chains(ref_pdb, subtype)
    ref_chain = cfg["reference_chain"]
    if ref_chain is None:
        if not ref_lbd:
            log(f"[FATAL] no PPAR LBD chain found in reference {ref_pdb}")
            _log_file.close()
            return
        ref_chain = ref_lbd[0]
    elif ref_chain not in ref_lbd:
        log(f"[WARN] reference chain {ref_chain} is not a PPAR LBD chain "
            f"by {method} check (LBD chains: {ref_lbd})")
    REF_OBJ = "reference"
    cmd.create(REF_OBJ, chain_sel(ref_pdb, ref_chain))
    cmd.delete(ref_pdb)
    ref_waters_sel = "ref_waters"
    cmd.select(ref_waters_sel, f"{REF_OBJ} and resn HOH")
    log(f"Reference {ref_pdb} chain {ref_chain} ({n_ca(REF_OBJ, ref_chain)} res, "
        f"{cmd.count_atoms(ref_waters_sel)} waters)\n")

    rows, chain_rows = [], []
    for pdb_id in cfg["pdb_ids"]:
        pdb_id = pdb_id.lower()
        if pdb_id == ref_pdb:
            log(f"Skipping {pdb_id}: reference structure")
            continue
        log(f"Processing {pdb_id}")
        try:
            load_structure(pdb_id, subtype)
        except Exception as exc:
            log(f"  [ERROR] could not load {pdb_id}: {exc}")
            chain_rows.append({"pdb_id": pdb_id, "all_protein_chains": "",
                               "lbd_chains": "", "used_chains": "",
                               "method": "load_failed", "note": str(exc)})
            continue

        all_chains = [c for c in cmd.get_chains(pdb_id)
                      if cmd.count_atoms(f"{pdb_id} and chain {c} and polymer.protein") > 0]
        lbd, method, note = identify_lbd_chains(pdb_id, subtype)
        if not lbd:
            log(f"  [SKIP] no PPAR LBD chain ({note})")
        else:
            if note:
                log(f"  [{method}] {note}")
        used = pick_chains(pdb_id, lbd, REF_OBJ) if lbd else []
        log(f"  protein chains {all_chains} -> LBD {lbd} -> using {used}")

        for c in used:
            align_and_collect(pdb_id, c, REF_OBJ, ref_waters_sel, rows)

        chain_rows.append({
            "pdb_id": pdb_id,
            "all_protein_chains": " ".join(all_chains),
            "lbd_chains": " ".join(lbd),
            "used_chains": " ".join(used),
            "method": method,
            "note": note,
        })
        cmd.delete(pdb_id)

    fieldnames = ["pdb_id", "chain", "water_resi", "x", "y", "z",
                  "align_rmsd", "n_aligned_atoms",
                  "align_rmsd_before_refine", "n_atoms_before_refine"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    with open(chains_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "all_protein_chains", "lbd_chains",
                                          "used_chains", "method", "note"])
        w.writeheader()
        w.writerows(chain_rows)

    n_len = sum(r["method"] == "length" for r in chain_rows)
    n_none = sum(r["method"] in ("none", "load_failed") for r in chain_rows)
    log(f"\nWrote {len(rows)} matched-water rows to {csv_path}")
    log(f"Wrote chain report to {chains_csv}")
    log(f"Entries using length fallback (check by hand): {n_len}")
    log(f"Entries skipped (no LBD / load failed): {n_none}")
    _log_file.close()
    _log_file = None


def main():
    for subtype in SUBTYPES_TO_RUN:
        run_subtype(subtype)


main()