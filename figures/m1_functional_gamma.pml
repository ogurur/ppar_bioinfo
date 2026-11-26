# Conserved Method 1 clusters (>=30%) and functional regions on 1PRG
# run from the project folder:  @figures/m1_functional_gamma.pml
reinitialize
load data/raw/PPARgamma/1prg.cif, ref
remove ref and not chain A
remove ref and resn HOH
bg_color white
hide everything
show cartoon, ref
color grey90, ref
set cartoon_transparency, 0.2
select pocket, ref and resi 288+289+323+367+449+473
select af2, ref and resi 293+294+297+298+301+304+306+311+312+314+315+318+319+466+467+468+469+470+471+472+473+474+475+476
select rxr, ref and resi 373+395+396+398+403+407+410+411+414+418+429+430+432+433+434+436+437+439+440+441+443+444+447+451
select selectivity, ref and resi 253+354
color lightorange, af2
color palecyan, rxr
show sticks, pocket and not name N+C+O
color orange, pocket and elem C
show sticks, selectivity and not name N+C+O
color violet, selectivity and elem C
load figures/method1_clusters_gamma.pdb, sites
hide everything, sites
set sphere_scale, 0.4
show spheres, sites and resi 1
color red, sites and resi 1
label sites and resi 1, "C1"
show spheres, sites and resi 2
color red, sites and resi 2
show spheres, sites and resi 3
color red, sites and resi 3
show spheres, sites and resi 4
color red, sites and resi 4
show spheres, sites and resi 5
color red, sites and resi 5
show spheres, sites and resi 6
color red, sites and resi 6
show spheres, sites and resi 7
color red, sites and resi 7
show spheres, sites and resi 8
color red, sites and resi 8
show spheres, sites and resi 9
color red, sites and resi 9
label sites and resi 9, "C9"
show spheres, sites and resi 10
color yelloworange, sites and resi 10
label sites and resi 10, "C10"
show spheres, sites and resi 11
color yelloworange, sites and resi 11
show spheres, sites and resi 12
color yelloworange, sites and resi 12
show spheres, sites and resi 13
color yelloworange, sites and resi 13
label sites and resi 13, "C13"
show spheres, sites and resi 14
color yelloworange, sites and resi 14
show spheres, sites and resi 15
color yelloworange, sites and resi 15
show spheres, sites and resi 16
color yelloworange, sites and resi 16
label sites and resi 16, "C16"
show spheres, sites and resi 17
color yelloworange, sites and resi 17
label sites and resi 17, "C17"
show spheres, sites and resi 18
color yelloworange, sites and resi 18
label sites and resi 18, "C18"
show spheres, sites and resi 19
color yelloworange, sites and resi 19
show spheres, sites and resi 20
color yelloworange, sites and resi 20
show spheres, sites and resi 21
color yelloworange, sites and resi 21
show spheres, sites and resi 22
color yelloworange, sites and resi 22
label sites and resi 22, "C22"
show spheres, sites and resi 23
color yelloworange, sites and resi 23
show spheres, sites and resi 24
color yelloworange, sites and resi 24
show spheres, sites and resi 25
color yelloworange, sites and resi 25
label sites and resi 25, "C25"
show spheres, sites and resi 26
color yelloworange, sites and resi 26
label sites and resi 26, "C26"
show spheres, sites and resi 27
color yelloworange, sites and resi 27
label sites and resi 27, "C27"
show spheres, sites and resi 28
color yelloworange, sites and resi 28
label sites and resi 28, "C28"
show spheres, sites and resi 29
color yelloworange, sites and resi 29
label sites and resi 29, "C29"
show spheres, sites and resi 30
color yelloworange, sites and resi 30
show spheres, sites and resi 31
color yelloworange, sites and resi 31
label sites and resi 31, "C31"
show spheres, sites and resi 32
color yelloworange, sites and resi 32
show spheres, sites and resi 33
color yelloworange, sites and resi 33
label sites and resi 33, "C33"
show spheres, sites and resi 34
color yelloworange, sites and resi 34
show spheres, sites and resi 35
color yelloworange, sites and resi 35
show spheres, sites and resi 36
color yelloworange, sites and resi 36
show spheres, sites and resi 38
color yelloworange, sites and resi 38
label sites and resi 38, "C38"
show spheres, sites and resi 37
color yelloworange, sites and resi 37
label sites and resi 37, "C37"
show spheres, sites and resi 40
color yelloworange, sites and resi 40
label sites and resi 40, "C40"
show spheres, sites and resi 39
color yelloworange, sites and resi 39
label sites and resi 39, "C39"
show spheres, sites and resi 42
color yelloworange, sites and resi 42
show spheres, sites and resi 41
color yelloworange, sites and resi 41
label sites and resi 41, "C41"
show spheres, sites and resi 43
color yelloworange, sites and resi 43
show spheres, sites and resi 44
color yelloworange, sites and resi 44
label sites and resi 44, "C44"
show spheres, sites and resi 45
color yelloworange, sites and resi 45
show spheres, sites and resi 46
color yelloworange, sites and resi 46
show spheres, sites and resi 47
color yelloworange, sites and resi 47
show spheres, sites and resi 48
color yelloworange, sites and resi 48
label sites and resi 48, "C48"
show spheres, sites and resi 49
color yelloworange, sites and resi 49
show spheres, sites and resi 50
color yelloworange, sites and resi 50
show spheres, sites and resi 51
color yelloworange, sites and resi 51
show spheres, sites and resi 52
color yelloworange, sites and resi 52
show spheres, sites and resi 53
color yelloworange, sites and resi 53
show spheres, sites and resi 54
color yelloworange, sites and resi 54
label sites and resi 54, "C54"
show spheres, sites and resi 55
color yelloworange, sites and resi 55
show spheres, sites and resi 56
color yelloworange, sites and resi 56
set label_color, black
deselect
orient ref
# orange sticks = pocket, light orange = AF-2 / coactivator groove,
# pale cyan = RXR interface, violet = selectivity residues;
# red spheres = high (>70%), yellow-orange = moderate (30-70%), labelled if near a region
