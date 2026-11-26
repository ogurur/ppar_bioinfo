# =====================================================================
# PPARgamma (1PRG chain A): top 5 conserved waters, the helices they are
# attached to, and their contact distances.
#  - helix names and colours as in the PPARalpha script
#  - boundaries = supervisor's table, PPARgamma1 numbering
#  - contacts = protein N/O atoms within 3.5 A of the water O (dashed, in A)
# Run in PyMOL from the project folder:  @figures/water_helices_gamma.pml
# Scenes: F1 overview, F2-F6 close-ups of W1-W5 (or: scene W1, recall)
# =====================================================================
reinitialize
fetch 1prg, async=0
create ppar_gamma, 1prg and chain A
delete 1prg
hide everything
show cartoon, ppar_gamma
color gray70, ppar_gamma
bg_color white
set cartoon_transparency, 0.15
set label_size, -0.8
set label_font_id, 7
set label_color, black
set float_labels, on
set label_distance_digits, 2
set dash_color, black
set dash_gap, 0.25
set dash_radius, 0.06
set ray_opaque_background, 0

# --- secondary-structure elements ---
select H1, ppar_gamma and polymer and resi 208-226
select H2, ppar_gamma and polymer and resi 230-238
select B1, ppar_gamma and polymer and resi 246-250
select H2p, ppar_gamma and polymer and resi 252-261
select H3, ppar_gamma and polymer and resi 274-300
select H3p, ppar_gamma and polymer and resi 306-309
select H4, ppar_gamma and polymer and resi 311-320
select H5, ppar_gamma and polymer and resi 322-333
select B2, ppar_gamma and polymer and resi 334-338
select B3, ppar_gamma and polymer and resi 339-342
select B4, ppar_gamma and polymer and resi 346-349
select H6, ppar_gamma and polymer and resi 350-355
select H7, ppar_gamma and polymer and resi 361-378
select H8, ppar_gamma and polymer and resi 380-394
select H9, ppar_gamma and polymer and resi 403-424
select H10, ppar_gamma and polymer and resi 430-457
select H12, ppar_gamma and polymer and resi 467-475
color tv_red, H1 and elem C
color smudge, H2 and elem C
color tv_blue, B1 and elem C
color tv_yellow, H2p and elem C
color magenta, H3 and elem C
color cyan, H3p and elem C
color orange, H4 and elem C
color pink, H5 and elem C
color lime, B2 and elem C
color salmon, B3 and elem C
color teal, B4 and elem C
color wheat, H6 and elem C
color firebrick, H7 and elem C
color skyblue, H8 and elem C
color chocolate, H9 and elem C
color palegreen, H10 and elem C
color deepblue, H12 and elem C

# helix labels (pseudoatoms named lab_*, so the selections above are kept)
pseudoatom lab_H1, selection=H1 and name CA
label lab_H1, "H1"
pseudoatom lab_H2, selection=H2 and name CA
label lab_H2, "H2"
pseudoatom lab_B1, selection=B1 and name CA
label lab_B1, "B1"
pseudoatom lab_H2p, selection=H2p and name CA
label lab_H2p, "H2p"
pseudoatom lab_H3, selection=H3 and name CA
label lab_H3, "H3"
pseudoatom lab_H3p, selection=H3p and name CA
label lab_H3p, "H3p"
pseudoatom lab_H4, selection=H4 and name CA
label lab_H4, "H4"
pseudoatom lab_H5, selection=H5 and name CA
label lab_H5, "H5"
pseudoatom lab_B2, selection=B2 and name CA
label lab_B2, "B2"
pseudoatom lab_B3, selection=B3 and name CA
label lab_B3, "B3"
pseudoatom lab_B4, selection=B4 and name CA
label lab_B4, "B4"
pseudoatom lab_H6, selection=H6 and name CA
label lab_H6, "H6"
pseudoatom lab_H7, selection=H7 and name CA
label lab_H7, "H7"
pseudoatom lab_H8, selection=H8 and name CA
label lab_H8, "H8"
pseudoatom lab_H9, selection=H9 and name CA
label lab_H9, "H9"
pseudoatom lab_H10, selection=H10 and name CA
label lab_H10, "H10"
pseudoatom lab_H12, selection=H12 and name CA
label lab_H12, "H12"
group helix_labels, lab_H* lab_B*
set label_size, -2.2, lab_H1
set label_size, -2.2, lab_H2
set label_size, -2.2, lab_B1
set label_size, -2.2, lab_H2p
set label_size, -2.2, lab_H3
set label_size, -2.2, lab_H3p
set label_size, -2.2, lab_H4
set label_size, -2.2, lab_H5
set label_size, -2.2, lab_B2
set label_size, -2.2, lab_B3
set label_size, -2.2, lab_B4
set label_size, -2.2, lab_H6
set label_size, -2.2, lab_H7
set label_size, -2.2, lab_H8
set label_size, -2.2, lab_H9
set label_size, -2.2, lab_H10
set label_size, -2.2, lab_H12

# --- conserved waters ---
select W1, ppar_gamma and resn HOH and resi 489
select W2, ppar_gamma and resn HOH and resi 614
select W3, ppar_gamma and resn HOH and resi 537
select W4, ppar_gamma and resn HOH and resi 498
select W5, ppar_gamma and resn HOH and resi 479
select top_waters, W1 or W2 or W3 or W4 or W5
show spheres, top_waters
set sphere_scale, 0.5, top_waters
color aquamarine, top_waters
label W1, "W1"
label W2, "W2"
label W3, "W3"
label W4, "W4"
label W5, "W5"
set label_color, teal, top_waters

# big water labels for the overview (pseudoatoms, like the helix labels)
pseudoatom lab_W1, selection=W1
label lab_W1, "W1"
set label_size, -2.2, lab_W1
set label_color, deepteal, lab_W1
set label_position, [0, 0, 3], lab_W1
pseudoatom lab_W2, selection=W2
label lab_W2, "W2"
set label_size, -2.2, lab_W2
set label_color, deepteal, lab_W2
set label_position, [0, 0, 3], lab_W2
pseudoatom lab_W3, selection=W3
label lab_W3, "W3"
set label_size, -2.2, lab_W3
set label_color, deepteal, lab_W3
set label_position, [0, 0, 3], lab_W3
pseudoatom lab_W4, selection=W4
label lab_W4, "W4"
set label_size, -2.2, lab_W4
set label_color, deepteal, lab_W4
set label_position, [0, 0, 3], lab_W4
pseudoatom lab_W5, selection=W5
label lab_W5, "W5"
set label_size, -2.2, lab_W5
set label_color, deepteal, lab_W5
set label_position, [0, 0, 3], lab_W5
group water_labels, lab_W*

# --- residues holding each water, labelled with residue and element ---
select c_W1, ppar_gamma and polymer and resi 301+303+306+307
label ppar_gamma and resi 301 and name CA, "Lys301 (loop)"
label ppar_gamma and resi 303 and name CA, "Ile303 (loop)"
label ppar_gamma and resi 306 and name CA, "Phe306 (H3p)"
label ppar_gamma and resi 307 and name CA, "Val307 (H3p)"
select c_W2, ppar_gamma and polymer and resi 381
label ppar_gamma and resi 381 and name CA, "Asp381 (H8)"
select c_W3, ppar_gamma and polymer and resi 230+375
label ppar_gamma and resi 230 and name CA, "Lys230 (H2)"
label ppar_gamma and resi 375 and name CA, "Asn375 (H7)"
select c_W4, ppar_gamma and polymer and resi 222+332
label ppar_gamma and resi 222 and name CA, "Tyr222 (H1)"
label ppar_gamma and resi 332 and name CA, "Ser332 (H5)"
select c_W5, ppar_gamma and polymer and resi 380+382
label ppar_gamma and resi 380 and name CA, "Asp380 (H8)"
label ppar_gamma and resi 382 and name CA, "Ser382 (H8)"
select contact_res, c_W1 or c_W2 or c_W3 or c_W4 or c_W5
util.cnc contact_res

# --- distances: water O to protein N/O within 3.5 A ---
# ILE303 O [loop] 2.47 A
distance d_W1_1, W1, ppar_gamma and resi 303 and name O
# LYS301 O [loop] 2.67 A
distance d_W1_2, W1, ppar_gamma and resi 301 and name O
# VAL307 N [H3p] 3.0 A
distance d_W1_3, W1, ppar_gamma and resi 307 and name N
# PHE306 N [H3p] 3.42 A
distance d_W1_4, W1, ppar_gamma and resi 306 and name N
# ASP381 OD2 [H8] 3.09 A
distance d_W2_1, W2, ppar_gamma and resi 381 and name OD2
# ASP381 N [H8] 3.15 A
distance d_W2_2, W2, ppar_gamma and resi 381 and name N
# ASN375 O [H7] 2.49 A
distance d_W3_1, W3, ppar_gamma and resi 375 and name O
# LYS230 NZ [H2] 3.24 A
distance d_W3_2, W3, ppar_gamma and resi 230 and name NZ
# TYR222 OH [H1] 2.7 A
distance d_W4_1, W4, ppar_gamma and resi 222 and name OH
# SER332 OG [H5] 2.74 A
distance d_W4_2, W4, ppar_gamma and resi 332 and name OG
# ASP380 OD1 [H8] 3.09 A
distance d_W5_1, W5, ppar_gamma and resi 380 and name OD1
# SER382 N [H8] 3.15 A
distance d_W5_2, W5, ppar_gamma and resi 382 and name N
# SER382 OG [H8] 3.47 A
distance d_W5_3, W5, ppar_gamma and resi 382 and name OG

# water-water bridges (second-shell waters = small grey spheres, grey dashes)
# W2 - HOH627 2.6 A
distance dw_W2_627, W2, ppar_gamma and resn HOH and resi 627
# W2 - HOH594 3.01 A
distance dw_W2_594, W2, ppar_gamma and resn HOH and resi 594
# W3 - HOH600 2.52 A
distance dw_W3_600, W3, ppar_gamma and resn HOH and resi 600
# W3 - HOH655 3.42 A
distance dw_W3_655, W3, ppar_gamma and resn HOH and resi 655
# W5 - HOH594 2.78 A
distance dw_W5_594, W5, ppar_gamma and resn HOH and resi 594
# W5 - HOH648 3.46 A
distance dw_W5_648, W5, ppar_gamma and resn HOH and resi 648
select second_shell, ppar_gamma and resn HOH and resi 594+600+627+648+655
set sphere_scale, 0.25, second_shell
color gray50, second_shell
color gray50, dw_*
group distances, d_* dw_*

# --- scene 1: overview (all waters, helix labels, no distances) ---
hide labels, polymer
hide sticks
hide spheres, second_shell
enable distances
disable d_* dw_*
enable helix_labels
enable water_labels
show cartoon, ppar_gamma
orient ppar_gamma
scene F1, store, overview: helices and the 5 conserved waters

# --- scene W1: close-up ---
hide labels, polymer
hide sticks
hide spheres, second_shell
disable d_* dw_*
disable helix_labels
disable water_labels
hide cartoon
show cartoon, byres (ppar_gamma and polymer within 10 of W1)
show sticks, byres c_W1
show labels, c_W1 and name CA
enable d_W1_*
orient W1 or c_W1
zoom W1 or c_W1, 4
scene F2, store, W1
scene W1, store

# --- scene W2: close-up ---
hide labels, polymer
hide sticks
hide spheres, second_shell
disable d_* dw_*
disable helix_labels
disable water_labels
hide cartoon
show cartoon, byres (ppar_gamma and polymer within 10 of W2)
show sticks, byres c_W2
show labels, c_W2 and name CA
enable d_W2_*
enable dw_W2_*
show spheres, ppar_gamma and resn HOH and resi 627+594
orient W2 or c_W2
zoom W2 or c_W2, 4
scene F3, store, W2
scene W2, store

# --- scene W3: close-up ---
hide labels, polymer
hide sticks
hide spheres, second_shell
disable d_* dw_*
disable helix_labels
disable water_labels
hide cartoon
show cartoon, byres (ppar_gamma and polymer within 10 of W3)
show sticks, byres c_W3
show labels, c_W3 and name CA
enable d_W3_*
enable dw_W3_*
show spheres, ppar_gamma and resn HOH and resi 600+655
orient W3 or c_W3
zoom W3 or c_W3, 4
scene F4, store, W3
scene W3, store

# --- scene W4: close-up ---
hide labels, polymer
hide sticks
hide spheres, second_shell
disable d_* dw_*
disable helix_labels
disable water_labels
hide cartoon
show cartoon, byres (ppar_gamma and polymer within 10 of W4)
show sticks, byres c_W4
show labels, c_W4 and name CA
enable d_W4_*
orient W4 or c_W4
zoom W4 or c_W4, 4
scene F5, store, W4
scene W4, store

# --- scene W5: close-up ---
hide labels, polymer
hide sticks
hide spheres, second_shell
disable d_* dw_*
disable helix_labels
disable water_labels
hide cartoon
show cartoon, byres (ppar_gamma and polymer within 10 of W5)
show sticks, byres c_W5
show labels, c_W5 and name CA
enable d_W5_*
enable dw_W5_*
show spheres, ppar_gamma and resn HOH and resi 594+648
orient W5 or c_W5
zoom W5 or c_W5, 4
scene F6, store, W5
scene W5, store

scene F1, recall

# To save pictures (300 dpi):
#   scene F1, recall
#   png figures/water_helices_gamma_overview.png, width=2400, height=1800, dpi=300, ray=1
#   scene W1, recall
#   png figures/water_helices_gamma_W1.png, width=2400, height=1800, dpi=300, ray=1
#   scene W2, recall
#   png figures/water_helices_gamma_W2.png, width=2400, height=1800, dpi=300, ray=1
#   scene W3, recall
#   png figures/water_helices_gamma_W3.png, width=2400, height=1800, dpi=300, ray=1
#   scene W4, recall
#   png figures/water_helices_gamma_W4.png, width=2400, height=1800, dpi=300, ray=1
#   scene W5, recall
#   png figures/water_helices_gamma_W5.png, width=2400, height=1800, dpi=300, ray=1
