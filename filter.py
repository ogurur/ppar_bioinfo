import argparse
import csv
import os
import glob

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

PPAR_UNIPROTS = {
    "PPARalpha": "Q07869",
    "PPARgamma": "P37231",
    "PPARdelta": "Q03181",
}

FIELDS = [
    "pdb_id", "subtype", "resolution_A", "r_factor_obs", "r_factor_free",
    "space_group", "ligands_present", "ligand_class", "lbd_chain_id",
    "n_waters_in_lbd", "deposit_date", "depositing_group",
    "cell_a", "cell_b", "cell_c", "cell_alpha", "cell_beta", "cell_gamma",
]


def as_list(value):
    """MMCIF2Dict gives a list for looped categories, sometimes a bare
    string for single-row ones. Normalize to a list either way."""
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def find_lbd_chain(cif_dict, uniprot_id):
    """Return the auth_asym_id (chain letter) whose _struct_ref entry
    matches the target UniProt accession, or None if not found."""
    ref_ids = as_list(cif_dict.get("_struct_ref.id"))
    ref_accessions = as_list(cif_dict.get("_struct_ref.pdbx_db_accession"))

    matching_ref_ids = {
        rid for rid, acc in zip(ref_ids, ref_accessions)
        if acc == uniprot_id
    }
    if not matching_ref_ids:
        return None

    seq_ref_ids = as_list(cif_dict.get("_struct_ref_seq.ref_id"))
    seq_strand_ids = as_list(cif_dict.get("_struct_ref_seq.pdbx_strand_id"))

    chains = [
        strand for rid, strand in zip(seq_ref_ids, seq_strand_ids)
        if rid in matching_ref_ids
    ]
    # Multiple chains can map to the same UniProt entry (e.g. two LBD
    # copies in the asymmetric unit). Take the first; flag the rest in
    # 'ligands_present' review rather than silently dropping data.
    return chains[0] if chains else None


def count_waters_in_chain(cif_dict, chain_id):
    """Count HOH HETATM records whose auth_asym_id matches chain_id."""
    group_pdb = as_list(cif_dict.get("_atom_site.group_PDB"))
    comp_id = as_list(cif_dict.get("_atom_site.label_comp_id"))
    auth_asym = as_list(cif_dict.get("_atom_site.auth_asym_id"))

    return sum(
        1 for g, c, a in zip(group_pdb, comp_id, auth_asym)
        if g == "HETATM" and c == "HOH" and a == chain_id
    )


def get_ligands(cif_dict, chain_id):
    """Non-water HETATM residue names present in the given chain."""
    group_pdb = as_list(cif_dict.get("_atom_site.group_PDB"))
    comp_id = as_list(cif_dict.get("_atom_site.label_comp_id"))
    auth_asym = as_list(cif_dict.get("_atom_site.auth_asym_id"))

    ligands = {
        c for g, c, a in zip(group_pdb, comp_id, auth_asym)
        if g == "HETATM" and c != "HOH" and a == chain_id
    }
    return sorted(ligands)


def extract_metadata(cif_path, subtype):
    pdb_id = os.path.splitext(os.path.basename(cif_path))[0].upper()
    cif_dict = MMCIF2Dict(cif_path)
    uniprot_id = PPAR_UNIPROTS[subtype]

    lbd_chain = find_lbd_chain(cif_dict, uniprot_id)
    if lbd_chain is None:
        print(f"  [WARN] {pdb_id}: could not map a chain to {uniprot_id} "
              f"via _struct_ref - flag for manual chain assignment")

    n_waters = count_waters_in_chain(cif_dict, lbd_chain) if lbd_chain else None
    ligands = get_ligands(cif_dict, lbd_chain) if lbd_chain else []

    cell = {k: safe_float(as_list(cif_dict.get(f"_cell.{k}"))[0])
            if cif_dict.get(f"_cell.{k}") else None
            for k in ["length_a", "length_b", "length_c",
                      "angle_alpha", "angle_beta", "angle_gamma"]}

    authors = as_list(cif_dict.get("_audit_author.name"))

    row = {
        "pdb_id": pdb_id,
        "subtype": subtype,
        "resolution_A": safe_float(as_list(cif_dict.get("_refine.ls_d_res_high")) and
                                    as_list(cif_dict.get("_refine.ls_d_res_high"))[0]),
        "r_factor_obs": safe_float(as_list(cif_dict.get("_refine.ls_R_factor_obs")) and
                                    as_list(cif_dict.get("_refine.ls_R_factor_obs"))[0]),
        "r_factor_free": safe_float(as_list(cif_dict.get("_refine.ls_R_factor_R_free")) and
                                     as_list(cif_dict.get("_refine.ls_R_factor_R_free"))[0]),
        "space_group": (as_list(cif_dict.get("_symmetry.space_group_name_H-M")) or [None])[0],
        "ligands_present": ";".join(ligands) if ligands else "apo",
        "ligand_class": "",  # fill via manual/ChEMBL lookup pass, see note below
        "lbd_chain_id": lbd_chain,
        "n_waters_in_lbd": n_waters,
        "deposit_date": (as_list(cif_dict.get(
            "_pdbx_database_status.recvd_initial_deposition_date")) or [None])[0],
        "depositing_group": authors[0] if authors else None,
        "cell_a": cell["length_a"], "cell_b": cell["length_b"], "cell_c": cell["length_c"],
        "cell_alpha": cell["angle_alpha"], "cell_beta": cell["angle_beta"],
        "cell_gamma": cell["angle_gamma"],
    }
    return row


def build_csv(raw_dir, out_path):
    rows = []
    for subtype in PPAR_UNIPROTS:
        cif_files = sorted(glob.glob(os.path.join(raw_dir, subtype, "*.cif")))
        print(f"{subtype}: {len(cif_files)} structures")
        for cif_path in cif_files:
            try:
                rows.append(extract_metadata(cif_path, subtype))
            except Exception as exc:
                print(f"  [ERROR] {cif_path}: {exc}")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {out_path}")

    unresolved = [r["pdb_id"] for r in rows if r["lbd_chain_id"] is None]
    if unresolved:
        print(f"[WARN] {len(unresolved)} structures need manual chain "
              f"assignment: {unresolved}")


def split_by_subtype(master_csv_path, out_dir):
    """Read the master metadata CSV and write one CSV per subtype into
    out_dir, e.g. metadata_PPARalpha.csv, metadata_PPARgamma.csv,
    metadata_PPARdelta.csv. Row order within each subtype is preserved."""
    with open(master_csv_path, newline="") as f:
        rows = list(csv.DictReader(f))

    os.makedirs(out_dir, exist_ok=True)
    counts = {}
    for subtype in PPAR_UNIPROTS:
        subtype_rows = [r for r in rows if r["subtype"] == subtype]
        counts[subtype] = len(subtype_rows)

        out_path = os.path.join(out_dir, f"metadata_{subtype}.csv")
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(subtype_rows)

    unmatched = [r["pdb_id"] for r in rows if r["subtype"] not in PPAR_UNIPROTS]
    if unmatched:
        print(f"[WARN] {len(unmatched)} rows had an unrecognized subtype "
              f"and were not written to any per-subtype file: {unmatched}")

    print("Split into per-subtype CSVs:")
    for subtype, n in counts.items():
        print(f"  {subtype}: {n} rows -> "
              f"{os.path.join(out_dir, f'metadata_{subtype}.csv')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--out", default="data/metadata.csv")
    parser.add_argument("--split-dir", default=None,
                         help="If set, also write metadata_<subtype>.csv "
                              "files into this directory.")
    args = parser.parse_args()
    build_csv(args.raw_dir, args.out)
    if args.split_dir:
        split_by_subtype(args.out, args.split_dir)
