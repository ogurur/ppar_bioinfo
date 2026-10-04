# Example: query PDB for PPARgamma structures with crystallographic waters
# Use the RCSB Search API (https://search.rcsb.org)

import requests, json
import os
from Bio.PDB import MMCIFParser, PDBList

def query_pdb_ppar(uniprot_id, min_resolution=None):
    """
    Query RCSB PDB for all X-ray structures of a given UniProt entry
    that contain at least one HOH HETATM record.
    uniprot_id: e.g. "P37231" (PPARgamma), "Q07869" (PPARalpha), "Q03181" (PPARdelta)
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                { "type": "terminal", "service": "text",
                  "parameters": { "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                                  "operator": "exact_match", "value": uniprot_id } },
                { "type": "terminal", "service": "text",
                  "parameters": { "attribute": "exptl.method",
                                  "operator": "exact_match", "value": "X-RAY DIFFRACTION" } },
            ]
        },
        "return_type": "entry",
        "request_options": { "return_all_hits": True }
    }
    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query",
                      json=query)
    return [hit["identifier"] for hit in r.json()["result_set"]]

def download_structures(pdb_ids, subtype, output_dir="data/raw"):
    """Download mmCIF files for a list of PDB IDs."""
    pdbl = PDBList()
    subdir = os.path.join(output_dir, subtype)
    os.makedirs(subdir, exist_ok=True)
    for pdb_id in pdb_ids:
        pdbl.retrieve_pdb_file(pdb_id.lower(), file_format="mmCif", pdir=subdir)
    print(f"Downloaded {len(pdb_ids)} structures for {subtype}")

ppar_uniprots = {
    "PPARalpha": "Q07869",
    "PPARgamma": "P37231",
    "PPARdelta": "Q03181",
}

all_ids = {}
for subtype, uid in ppar_uniprots.items():
    ids = query_pdb_ppar(uid)

    print(f"{subtype}: {len(ids)} structures found")
    all_ids[subtype] = ids

# print(all_ids)

#print(all_ids['PPARalpha'])
#download_structures(all_ids["PPARalpha"], 'PPARalpha')

# print(all_ids['PPARdelta'])
# download_structures(all_ids["PPARdelta"], 'PPARdelta')

print(all_ids['PPARgamma'])
#download_structures(all_ids["PPARgamma"], 'PPARgamma')
