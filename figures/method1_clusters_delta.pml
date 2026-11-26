# Method 1 clusters (>=30%) on the 1GWX reference, sphere colour and size = conservation
# run from the project folder:  @figures/method1_clusters_delta.pml
reinitialize
load data/raw/PPARdelta/1gwx.cif, ref
remove ref and not chain A
remove ref and resn HOH
hide everything
show cartoon, ref
color grey80, ref
bg_color white
load figures/method1_clusters_delta.pdb, sites
show spheres, sites
set sphere_scale, 0.35, sites
set sphere_scale, 0.55, sites and b > 70
spectrum b, yellow_red, sites, minimum=30, maximum=90
label sites and b > 70, "C%s" % resi.strip()
set label_color, black
orient ref
