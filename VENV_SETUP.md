# Virtual Environment Setup for PPAR Bioinfo Project

## Summary
A Python virtual environment has been created with all necessary dependencies for the project.

## Installed Packages
- **biopython** 1.88 - Bioinformatics toolkit for structural biology (PDB parsing, MMCIF handling)
- **numpy** 2.2.6 - Numerical computing
- **pandas** 2.3.3 - Data manipulation and CSV handling
- **scipy** 1.15.3 - Scientific computing (hierarchical clustering)
- **requests** 2.33.1 - HTTP library for API requests

## Project Dependencies Scanned
The following modules are imported in this project:
- `Bio.PDB` (BioPython) - Protein Data Bank file parsing
- `pymol` - ⚠️ Not available via pip (requires Schrodinger installation)
- `scipy.cluster.hierarchy` - Hierarchical clustering algorithms
- Standard library: argparse, csv, glob, json, os

## Usage Instructions

### Option 1: Use the created venv (Recommended)
The virtual environment has been created in `/sessions/beautiful-upbeat-knuth/mnt/outputs/ppar_venv`

```bash
# Activate the environment
source /sessions/beautiful-upbeat-knuth/mnt/outputs/ppar_venv/bin/activate

# Run your scripts
python alignment_pparalpha.py
```

### Option 2: Create a local venv in the project folder
```bash
cd ~/path/to/ppar_bioinfo
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Note on PyMOL
PyMOL is not available through standard pip. It requires installation from Schrodinger:
- Visit: https://pymol.org/
- Install following Schrodinger's instructions
- Once installed, the scripts should work without modification

## Files Generated
- `requirements.txt` - List of pip-installable dependencies
- `VENV_SETUP.md` - This file
