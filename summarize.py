import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
 
# subtype -> (csv file, number of structures processed, i.e. excluding reference)
SUBTYPES = {
    "alpha": ("conserved_waters_alpha.csv", 78),
    "gamma": ("conserved_waters_gamma.csv", 378),
    "delta": ("conserved_waters_delta.csv", 54),
}
SITE_CUTOFF = 1.0  # Angstroms
 
summary_rows = []
for subtype, (path, n_structures) in SUBTYPES.items():
    df = pd.read_csv(path, dtype={"pdb_id": str})
 
    per_struct = df.groupby("pdb_id").size()
    n_zero = n_structures - len(per_struct)
    counts = np.concatenate([per_struct.values, np.zeros(n_zero)])
 
    chains = df.drop_duplicates(["pdb_id", "chain"])
 
    labels = fcluster(linkage(df[["x", "y", "z"]].values, "average"),
                      SITE_CUTOFF, "distance")
    df["site"] = labels
    site_freq = df.groupby("site")["pdb_id"].nunique() / n_structures
    site_freq = site_freq.sort_values(ascending=False)
    site_freq.rename("fraction_of_structures").to_csv(f"water_sites_{subtype}.csv")
 
    summary_rows.append({
        "subtype": subtype,
        "structures_examined": n_structures,
        "matched_water_rows": len(df),
        "mean_per_structure": round(counts.mean(), 1),
        "sd_per_structure": round(counts.std(ddof=1), 1),
        "median_per_structure": float(np.median(counts)),
        "max_per_structure": f"{per_struct.max()} ({per_struct.idxmax()})",
        "structures_with_0": n_zero,
        "mean_align_rmsd": round(chains["align_rmsd"].mean(), 2),
        "n_sites": len(site_freq),
        "sites_ge_50pct": int((site_freq >= 0.5).sum()),
        "sites_ge_25pct": int((site_freq >= 0.25).sum()),
        "top_site_pct": round(100 * site_freq.iloc[0]),
    })
 
summary = pd.DataFrame(summary_rows)
print(summary.to_string(index=False))
summary.to_csv("water_summary.csv", index=False)