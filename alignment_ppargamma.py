# if you have .cif files on your computer - write directory in LOCAL_CIF_DIR like this "C:/Users/Ogurur/ppar/data/raw/PPARdelta" 
# if not - just keep as it is

import os
import csv
from pymol import cmd

# ---------------------------------------------------------------------
# Config - edit these for your run
# ---------------------------------------------------------------------
REFERENCE_PDB = "1prg" # 2prg
REFERENCE_CHAIN = "A"
REF_OBJ = "reference" 

CUTOFF = 1.0                     # Angstroms, same as your `within 1.0`
OUTPUT_CSV = "conserved_waters_gamma.csv"
LOG_TXT = "conserved_waters_log_gamma.txt"
LOCAL_CIF_DIR = "C:/Users/Ogurur/ppar/data/raw/PPARgamma"        # e.g. "data/raw/PPARgamma" - if set,
                                  # loads local files instead of cmd.fetch

# Structures to process. Add every PDB ID from your run here.
PDB_IDS = [
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
]


_log_file = None  # opened in main(), closed at the end

def log(message=""):
    """Print to the console AND append to LOG_TXT, so the run's output
    is saved even if you can't scroll back far enough in PyMOL."""
    print(message)
    if _log_file is not None:
        _log_file.write(str(message) + "\n")
        _log_file.flush()  # write immediately in case PyMOL crashes mid-run


def load_structure(pdb_id):
    """Load a structure either from a local mmCIF or via PyMOL's fetch."""
    if LOCAL_CIF_DIR:
        path = os.path.join(LOCAL_CIF_DIR, f"{pdb_id}.cif")
        cmd.load(path, pdb_id)
    else:
        cmd.fetch(pdb_id, type="cif")
 
 
def get_polymer_chains(obj_name):
    """Return chain letters that actually contain polymer (protein),
    skipping waters/ligand-only chains."""
    chains = cmd.get_chains(obj_name)
    protein_chains = [
        c for c in chains
        if cmd.count_atoms(f"{obj_name} and chain {c} and polymer.protein") > 0
    ]
    return protein_chains if protein_chains else [""]
 
 
def align_and_collect(obj_name, chain, reference, ref_waters_sel,
                       cutoff, rows, pdb_id):
    """Align obj_name (optionally restricted to one chain) onto the
    reference, then record any waters within `cutoff` A of the
    reference's waters."""
    work_obj = obj_name
    if chain:
        work_obj = f"{obj_name}_{chain}"
        cmd.create(work_obj, f"{obj_name} and chain {chain}")
 
    n_align = cmd.get_model(work_obj).nAtom
    if n_align == 0:
        log(f"  [skip] {work_obj}: no atoms")
        return
 
    try:
        align_result = cmd.align(work_obj, reference)
    except Exception as exc:
        log(f"  [ERROR] {work_obj} alignment failed: {exc}")
        return
 
    # cmd.align returns (RMSD after refinement, n atoms after refinement, ...)
    align_rmsd = align_result[0]
    n_aligned = align_result[1]
 
    sel_name = f"matched_{work_obj}"
    cmd.select(
        sel_name,
        f"({work_obj} and resn HOH) within {cutoff} of {ref_waters_sel}"
    )
    n_matched = cmd.count_atoms(sel_name)
 
    model = cmd.get_model(sel_name)
    for atom in model.atom:
        rows.append({
            "pdb_id": pdb_id,
            "chain": chain or "single",
            "water_resi": atom.resi,
            "x": round(atom.coord[0], 2),
            "y": round(atom.coord[1], 2),
            "z": round(atom.coord[2], 2),
            "align_rmsd": round(align_rmsd, 3),
            "n_aligned_atoms": n_aligned,
        })
 
    log(f"  {work_obj}: {n_matched} matched waters "
          f"(align RMSD={align_rmsd:.3f} over {n_aligned} atoms)")
 
    if chain:
        cmd.delete(work_obj)
    cmd.delete(sel_name)
 
 
def main():
    global _log_file
    log_path = os.path.abspath(LOG_TXT)
    csv_path = os.path.abspath(OUTPUT_CSV)
    _log_file = open(log_path, "w")
 
    cmd.reinitialize()
 
    log(f"Current working directory: {os.getcwd()}")
    log(f"Log will be saved to: {log_path}")
    log(f"CSV will be saved to: {csv_path}")
    log(f"Loading reference: {REFERENCE_PDB} chain {REFERENCE_CHAIN}")
    load_structure(REFERENCE_PDB)
    log(f"  chains in {REFERENCE_PDB}: {cmd.get_chains(REFERENCE_PDB)}")
    cmd.create(REF_OBJ, f"{REFERENCE_PDB} and chain {REFERENCE_CHAIN}")
    cmd.delete(REFERENCE_PDB)
    ref_waters_sel = "ref_waters"
    cmd.select(ref_waters_sel, f"{REF_OBJ} and resn HOH")
    n_ref_waters = cmd.count_atoms(ref_waters_sel)
    log(f"Reference (single chain) has {n_ref_waters} water atoms\n")
 
    rows = []
    for pdb_id in PDB_IDS:
        if pdb_id.lower() == REFERENCE_PDB.lower():
            log(f"Skipping {pdb_id}: it's the reference structure, "
                f"already loaded above")
            continue
 
        log(f"Processing {pdb_id}")
        try:
            load_structure(pdb_id)
        except Exception as exc:
            log(f"  [ERROR] could not load {pdb_id}: {exc}")
            continue
 
        chains = get_polymer_chains(pdb_id)
        # Single-chain structures: chains might be [''] or ['A'] with
        # only one real chain. If there's exactly one chain, align the
        # whole object directly, matching how you handled 1y0s / 2gwx etc.
        if len(chains) <= 1:
            align_and_collect(pdb_id, None, REF_OBJ,
                               ref_waters_sel, CUTOFF, rows, pdb_id)
        else:
            for chain in chains:
                align_and_collect(pdb_id, chain, REF_OBJ,
                                   ref_waters_sel, CUTOFF, rows, pdb_id)
 
        cmd.delete(pdb_id)
 
    fieldnames = ["pdb_id", "chain", "water_resi", "x", "y", "z",
                  "align_rmsd", "n_aligned_atoms"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
 
    log(f"\nWrote {len(rows)} matched-water rows to {csv_path}")
    log(f"Saved run log to {log_path}")
    _log_file.close()
 
 
main()