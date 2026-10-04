# if you have files on your computer - write directory in LOCAL_CIF_DIR like this "C:/Users/Ogurur/ppar/data/raw/PPARdelta" 
# if not - just keep as it is

import os
import csv
from pymol import cmd

# ---------------------------------------------------------------------
# Config - edit these for your run
# ---------------------------------------------------------------------
REFERENCE_PDB = "1k7l"
CUTOFF = 1.0                     # Angstroms, same as your `within 1.0`
OUTPUT_CSV = "conserved_waters_alpha.csv"
LOG_TXT = "conserved_waters_log_alpha.txt"
LOCAL_CIF_DIR = "C:/Users/Ogurur/ppar/data/raw/PPARalpha"          # e.g. "data/raw/PPARgamma" - if set,
                                  # loads local files instead of cmd.fetch

# Structures to process. Add every PDB ID from your run here.
PDB_IDS = [
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
    return chains if chains else [""]
 
 
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
        super_result = cmd.super(work_obj, reference)
    except Exception as exc:
        log(f"  [ERROR] {work_obj} alignment failed: {exc}")
        return
 
    align_rmsd = align_result[0]
    super_rmsd = super_result[0]
 
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
            "super_rmsd": round(super_rmsd, 3),
        })
 
    log(f"  {work_obj}: {n_matched} matched waters "
          f"(align RMSD={align_rmsd:.3f}, super RMSD={super_rmsd:.3f})")
 
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
    log(f"Loading reference: {REFERENCE_PDB}")
    load_structure(REFERENCE_PDB)
    ref_waters_sel = "ref_waters"
    cmd.select(ref_waters_sel, f"resn HOH and {REFERENCE_PDB}")
    n_ref_waters = cmd.count_atoms(ref_waters_sel)
    log(f"Reference has {n_ref_waters} water atoms\n")
 
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
            align_and_collect(pdb_id, None, REFERENCE_PDB,
                               ref_waters_sel, CUTOFF, rows, pdb_id)
        else:
            for chain in chains:
                align_and_collect(pdb_id, chain, REFERENCE_PDB,
                                   ref_waters_sel, CUTOFF, rows, pdb_id)
 
        cmd.delete(pdb_id)
 
    fieldnames = ["pdb_id", "chain", "water_resi", "x", "y", "z",
                  "align_rmsd", "super_rmsd"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
 
    log(f"\nWrote {len(rows)} matched-water rows to {csv_path}")
    log(f"Saved run log to {log_path}")
    _log_file.close()
 
 
main()
 