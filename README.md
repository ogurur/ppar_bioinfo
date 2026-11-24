# PPAR Water Conservation Analysis

BIOL363 Structural Bioinformatics Project 2026. Supervisor: F. Molnár, Nazarbayev University.
Team: Raniya Gaini, Nuradil Amangeldi, Nassiba Omarkhanova, Kamila Nurzhanova.

The project finds crystallographic water molecules that sit at the same position in most
PPARα, PPARγ and PPARδ ligand-binding-domain (LBD) structures in the PDB, and checks whether
they are real (ordered, resolution-independent, not one lab's repeats) and where they sit in
the protein. See `PPAR_Water_Conservation_Project (1).docx` for the full project description.

---

## Setup

Two environments are used:

| What | Environment | Install |
|---|---|---|
| All analysis scripts (everything except the PyMOL alignment) | `.venv` (or any Python 3.9+) | `pip install biopython numpy pandas scipy scikit-learn matplotlib requests` |
| `alignment_ppar_lbd*.py` and viewing `.pml` figures | PyMOL 3.x (conda `ppar` env or PyMOL GUI) | `conda install -c conda-forge pymol-open-source` |

Run every script **from the project folder** (`C:\Users\Ogurur\ppar_bioinfo`), e.g.
`python method1_prepare.py`. Scripts find their input files relative to that folder.
Inside PyMOL, run `.pml` files with `@figures/<name>.pml` (not `pymol <file>`).

---

## Folder layout

```
data/
  raw/PPARalpha|PPARgamma|PPARdelta/   downloaded mmCIF files, one per PDB entry
  metadata.csv                         resolution, R-factors, space group, ligands, cell, lab, date
  processed/method1/                   aligned-water and superposition tables (method1_prepare.py)
  results/                             result tables of the analysis scripts
  results/method1/                     Method 1 cluster atlas and sensitivity tables
figures/                               all plots (.png, 300 DPI) and PyMOL scripts (.pml, .pdb)
lbd_chains_PPAR<subtype>.csv           which chain of each entry is the PPAR LBD
conserved_align_waters_<subtype>[_05|_2].csv   reference-based matched waters (1 / 0.5 / 2 Å)
```

---

## Pipeline, in the order to run it

### 1. Dataset

| Script | What it does | Output |
|---|---|---|
| `script.py` | Queries the RCSB Search API for all X-ray structures of PPARα (Q07869), PPARγ (P37231) and PPARδ (Q03181); the download call (`download_structures`) is commented out — uncomment to download mmCIF files into `data/raw/<subtype>/`. | PDB ID lists, mmCIF files |
| `filter.py` | Reads every mmCIF header and builds the metadata table: resolution, R/R-free, space group, ligands, PPAR LBD chain (by UniProt accession), number of LBD waters, deposition date, depositing group, unit cell. Run `python filter.py` (options `--raw-dir`, `--out`; `--split-dir data/proccessed` also writes one file per subtype). | `data/metadata.csv` (+ `metadata_<subtype>.csv`) |

### 2. Reference-based water matching (PyMOL)

| Script | What it does | Output |
|---|---|---|
| `alignment_ppar_lbd.py` | Run **in PyMOL** (or with conda PyMOL's Python). For each subtype it keeps only the PPAR LBD chain (UniProt accession, length fallback), superimposes it on the reference (1K7L α, 1PRG γ, 1GWX δ) with `cmd.align`, and records every water within `CUTOFF` (1.0 Å) of a reference water. One chain per entry (lowest RMSD). | `conserved_align_waters_<subtype>.csv`, `lbd_chains_PPAR<subtype>.csv`, run logs |
| `alignment_ppar_lbd_05.py`, `alignment_ppar_lbd_2.py` | The same script with `CUTOFF` = 0.5 Å and 2.0 Å (cutoff sensitivity). | `..._05.csv`, `..._2.csv` |

### 3. Conservation of the reference-based sites

| Script | What it does | Output |
|---|---|---|
| `plot_water_conservation.py` | Assigns each matched water to its nearest reference water (a *site*), computes conservation = structures with a water at the site ÷ aligned structures, names sites W1, W2, … from most to least conserved, and plots them for 0.5 / 1 / 2 Å. | `data/results/water_sites_conservation_<subtype>.csv`, `conservation_summary.csv`, `figures/conservation_<subtype>_bars.png`, `..._ranked.png` |

### 4. Analyses of the conserved sites (supervisor's tasks)

| Script | Question | Output |
|---|---|---|
| `label_helices.py` | Task 1 – which helix / loop is each conserved site next to? Helices named in PPARγ1 numbering; α and δ mapped by sequence alignment. | `data/results/water_site_helices_<subtype>.csv`, `figures/sites_per_element.png`, `figures/pymol_sites_<subtype>.pml` |
| `water_bfactors.py` | Task 2 – are conserved waters more ordered? Water B-factor ÷ mean Cα B-factor of its chain, per site; Spearman and Mann–Whitney tests. | `data/results/water_site_bfactors_<subtype>.csv`, `bfactor_stats.csv`, `figures/bfactor_by_class.png`, `bfactor_vs_conservation.png` |
| `resolution_vs_conservation.py` | Task 3 – does conservation depend on resolution? Per structure and per site, with Spearman ρ. | `data/results/resolution_per_structure_<subtype>.csv`, `resolution_per_site_<subtype>.csv`, `resolution_stats.csv`, `figures/resolution_vs_conservation.png` |
| `plot_resolution_heatmap.py` | Heatmap of each site's occupancy in ≤2.0, 2.0–2.5, >2.5 Å structures (run after the script above). | `figures/resolution_heatmap.png` |
| `zero_match_structures.py` | Task 4 – which structures have no water at any site, and why (resolution, no modelled waters, lab)? | `data/results/zero_match_structures.csv`, `zero_match_stats.csv`, `figures/zero_match_resolution.png` |
| `water_contacts.py` | Contacts and bridges: residues each conserved water H-bonds (residue–H₂O–residue, residue–H₂O–H₂O–residue) in every structure, in reference numbering; cross-subtype comparison of sites. | `data/results/water_contacts_<subtype>.csv`, `water_site_contacts_<subtype>.csv`, `cross_subtype_sites.csv`, `cross_subtype_occupancy.csv` |
| `plot_cross_subtype.py` | Heatmap of which conserved sites are shared by α, γ and δ (run after `water_contacts.py`). | `figures/cross_subtype_heatmap.png` |

### 5. Method 1 – full version (project description, section 4)

Run in this order; no PyMOL needed.

| Script | What it does | Output |
|---|---|---|
| `method1_prepare.py` | Superimposes every LBD chain on its reference using **Cα of the core helices** (H3, H4/5, H7, H10/11; no outlier rejection) and also on the full LBD; does the same onto two alternative references per subtype; extracts **all** waters within 5 Å of the protein with B-factors. | `data/processed/method1/waters_<subtype>.csv`, `transforms_<subtype>.csv`, `references.csv` |
| `method1_cluster.py` | Clusters all waters (distance-cutoff clustering, main setting 1.5 Å, and DBSCAN); conservation per cluster, also counting near-duplicate structures (same lab, space group, unit cell) once; matches clusters to the reference-based sites; sensitivity to clustering cutoff (1.0–3.5 Å), resolution, alignment scope, RMSD outliers and reference choice. Large subtypes can be run in parts: `python method1_cluster.py . gamma sweep`, then `... gamma sens`. | `data/results/method1/clusters_<subtype>.csv` (conservation atlas), `near_duplicate_groups_<subtype>.csv`, `sensitivity_cutoff.csv`, `sensitivity_other.csv`; `figures/method1_rmsd.png`, `method1_sensitivity_cutoff.png`, `method1_clusters_ranked.png`, `method1_clusters_<subtype>.pdb/.pml` |

---

## Key settings

| Setting | Value | Where |
|---|---|---|
| References | 1K7L (α), 1PRG (γ), 1GWX (δ), chain A | all scripts (`SUBTYPES` / `CONFIG`) |
| Water cutoff (reference-based) | 1.0 Å main; 0.5 and 2.0 Å for sensitivity | `alignment_ppar_lbd*.py` |
| Alignment refinement (PyMOL `align`) | PyMOL default: 5 cycles, 2 Å | `alignment_ppar_lbd.py` |
| Core helices (Method 1) | PPARγ1 276–303, 306–333, 359–377, 431–459 | `method1_prepare.py` (`CORE`) |
| Clustering (Method 1) | distance cutoff 1.5 Å | `method1_cluster.py` (`MAIN`) |
| Conservation classes | high > 70 %, moderate 30–70 %, rare < 30 % | all analysis scripts |
| RMSD outlier | core Cα RMSD > 2.0 Å | `method1_cluster.py` (`MAX_RMSD`) |
| Contact / H-bond distance | 3.5 Å (contacts), 4.0 Å (helix proximity) | `water_contacts.py`, `label_helices.py` |

---

## Older scripts (kept for reference, superseded)

| Script | Status |
|---|---|
| `alignment_pparalpha.py`, `alignment_ppargamma.py`, `alignment_ppardelta.py` | First per-subtype PyMOL scripts (all chains, `cmd.super`). Replaced by `alignment_ppar_lbd.py`. Their outputs are `conserved_waters_<subtype>.csv`. |
| `summarize.py` | Summary of the old `conserved_waters_<subtype>.csv` files with hierarchical clustering at 1 Å → `water_sites_<subtype>.csv`, `water_summary.csv`. Replaced by `plot_water_conservation.py`. |
| `plot_step_a.py` | Early QC plots (RMSD histogram, waters per chain, density projections) from the old CSVs → `figures/01_…05_*.png`. |
| `test.py` | Checks that the raw-data folder path exists. |
| `VENV_SETUP.md` | Notes from the first virtual-environment setup. |

---

## Notes

- PyMOL and the Method 1 scripts align differently (PyMOL `align` on all LBD atoms vs Kabsch on
  core-helix Cα), so their percentages differ slightly; the top sites agree.
- Helix names H6–H9 were assigned by order in the chain; check them against a PPARγ structure
  paper before publishing.
- Structures with no modelled waters are excluded from Method 1 percentages; the
  reference-based percentages include them.
