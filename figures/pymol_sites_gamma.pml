# PPARγ: conserved water sites (>= 30% at 1 Å) and nearby helices
# run from the project folder:  pymol figures/pymol_sites_gamma.pml
reinitialize
load data/raw/PPARgamma/1prg.cif, 1prg
remove 1prg and not chain A
bg_color white
hide everything
show cartoon, 1prg
color grey90, 1prg
set cartoon_transparency, 0.15
set label_size, 16
set label_color, black
set float_labels, 1
select hx_H1, 1prg and resi 207+208+209+210+211+212+213+214+215+216+217+218+219+220+221+222+223+224+225
color paleyellow, hx_H1
pseudoatom lab_H1, hx_H1 and name CA
label lab_H1, "H1"
hide everything, lab_H1
show labels, lab_H1
select hx_H2, 1prg and resi 230+231+232+233+234+235+236+237
color grey70, hx_H2
pseudoatom lab_H2, hx_H2 and name CA
label lab_H2, "H2"
hide everything, lab_H2
show labels, lab_H2
select hx_H2p, 1prg and resi 252+253+254+255+256+257+258+259
color grey60, hx_H2p
pseudoatom lab_H2p, hx_H2p and name CA
label lab_H2p, "H2'"
hide everything, lab_H2p
show labels, lab_H2p
select hx_H3, 1prg and resi 276+277+278+279+280+281+282+283+284+285+286+287+288+289+290+291+292+293+294+295+296+297+298+299+300+301+302+303
color slate, hx_H3
pseudoatom lab_H3, hx_H3 and name CA
label lab_H3, "H3"
hide everything, lab_H3
show labels, lab_H3
select hx_H4_5, 1prg and resi 306+307+308+309+310+311+312+313+314+315+316+317+318+319+320+321+322+323+324+325+326+327+328+329+330+331+332+333
color aquamarine, hx_H4_5
pseudoatom lab_H4_5, hx_H4_5 and name CA
label lab_H4_5, "H4/5"
hide everything, lab_H4_5
show labels, lab_H4_5
select hx_H6, 1prg and resi 350+351+352+353+354+355
color wheat, hx_H6
pseudoatom lab_H6, hx_H6 and name CA
label lab_H6, "H6"
hide everything, lab_H6
show labels, lab_H6
select hx_H7, 1prg and resi 359+360+361+362+363+364+365+366+367+368+369+370+371+372+373+374+375+376+377
color lightorange, hx_H7
pseudoatom lab_H7, hx_H7 and name CA
label lab_H7, "H7"
hide everything, lab_H7
show labels, lab_H7
select hx_H8, 1prg and resi 381+382+383+384+385+386+387+388+389+390+391+392
color palegreen, hx_H8
pseudoatom lab_H8, hx_H8 and name CA
label lab_H8, "H8"
hide everything, lab_H8
show labels, lab_H8
select hx_H9, 1prg and resi 403+404+405+406+407+408+409+410+411+412+413+414+415+416+417+418+419+420+421+422+423+424
color lightpink, hx_H9
pseudoatom lab_H9, hx_H9 and name CA
label lab_H9, "H9"
hide everything, lab_H9
show labels, lab_H9
select hx_H10_11, 1prg and resi 431+432+433+434+435+436+437+438+439+440+441+442+443+444+445+446+447+448+449+450+451+452+453+454+455+456+457+458+459
color violet, hx_H10_11
pseudoatom lab_H10_11, hx_H10_11 and name CA
label lab_H10_11, "H10/11"
hide everything, lab_H10_11
show labels, lab_H10_11
select hx_H12, 1prg and resi 466+467+468+469+470+471+472+473+474+475+476
color salmon, hx_H12
pseudoatom lab_H12, hx_H12 and name CA
label lab_H12, "H12"
hide everything, lab_H12
show labels, lab_H12
select site_W1, 1prg and resn HOH and resi 489
show spheres, site_W1
set sphere_scale, 0.45, site_W1
color red, site_W1
label site_W1, "W1"
select site_W2, 1prg and resn HOH and resi 614
show spheres, site_W2
set sphere_scale, 0.45, site_W2
color red, site_W2
label site_W2, "W2"
select site_W3, 1prg and resn HOH and resi 537
show spheres, site_W3
set sphere_scale, 0.45, site_W3
color red, site_W3
label site_W3, "W3"
select site_W4, 1prg and resn HOH and resi 498
show spheres, site_W4
set sphere_scale, 0.45, site_W4
color red, site_W4
label site_W4, "W4"
select site_W5, 1prg and resn HOH and resi 479
show spheres, site_W5
set sphere_scale, 0.45, site_W5
color red, site_W5
label site_W5, "W5"
select site_W6, 1prg and resn HOH and resi 719
show spheres, site_W6
set sphere_scale, 0.3, site_W6
color orange, site_W6
label site_W6, "W6"
select site_W7, 1prg and resn HOH and resi 553
show spheres, site_W7
set sphere_scale, 0.3, site_W7
color orange, site_W7
label site_W7, "W7"
select site_W8, 1prg and resn HOH and resi 662
show spheres, site_W8
set sphere_scale, 0.3, site_W8
color orange, site_W8
label site_W8, "W8"
select site_W9, 1prg and resn HOH and resi 656
show spheres, site_W9
set sphere_scale, 0.3, site_W9
color orange, site_W9
label site_W9, "W9"
select site_W10, 1prg and resn HOH and resi 559
show spheres, site_W10
set sphere_scale, 0.3, site_W10
color orange, site_W10
label site_W10, "W10"
select site_W11, 1prg and resn HOH and resi 691
show spheres, site_W11
set sphere_scale, 0.3, site_W11
color orange, site_W11
label site_W11, "W11"
select site_W12, 1prg and resn HOH and resi 523
show spheres, site_W12
set sphere_scale, 0.3, site_W12
color orange, site_W12
label site_W12, "W12"
select site_W13, 1prg and resn HOH and resi 482
show spheres, site_W13
set sphere_scale, 0.3, site_W13
color orange, site_W13
label site_W13, "W13"
select site_W14, 1prg and resn HOH and resi 743
show spheres, site_W14
set sphere_scale, 0.3, site_W14
color orange, site_W14
label site_W14, "W14"
select site_W15, 1prg and resn HOH and resi 627
show spheres, site_W15
set sphere_scale, 0.3, site_W15
color orange, site_W15
label site_W15, "W15"
select site_W16, 1prg and resn HOH and resi 591
show spheres, site_W16
set sphere_scale, 0.3, site_W16
color orange, site_W16
label site_W16, "W16"
select site_W17, 1prg and resn HOH and resi 706
show spheres, site_W17
set sphere_scale, 0.3, site_W17
color orange, site_W17
label site_W17, "W17"
select site_W18, 1prg and resn HOH and resi 562
show spheres, site_W18
set sphere_scale, 0.3, site_W18
color orange, site_W18
label site_W18, "W18"
select site_W19, 1prg and resn HOH and resi 515
show spheres, site_W19
set sphere_scale, 0.3, site_W19
color orange, site_W19
label site_W19, "W19"
select site_W20, 1prg and resn HOH and resi 714
show spheres, site_W20
set sphere_scale, 0.3, site_W20
color orange, site_W20
label site_W20, "W20"
select site_W21, 1prg and resn HOH and resi 760
show spheres, site_W21
set sphere_scale, 0.3, site_W21
color orange, site_W21
label site_W21, "W21"
select site_W22, 1prg and resn HOH and resi 649
show spheres, site_W22
set sphere_scale, 0.3, site_W22
color orange, site_W22
label site_W22, "W22"
select site_W23, 1prg and resn HOH and resi 617
show spheres, site_W23
set sphere_scale, 0.3, site_W23
color orange, site_W23
label site_W23, "W23"
select site_W24, 1prg and resn HOH and resi 722
show spheres, site_W24
set sphere_scale, 0.3, site_W24
color orange, site_W24
label site_W24, "W24"
select site_W25, 1prg and resn HOH and resi 520
show spheres, site_W25
set sphere_scale, 0.3, site_W25
color orange, site_W25
label site_W25, "W25"
select site_W26, 1prg and resn HOH and resi 588
show spheres, site_W26
set sphere_scale, 0.3, site_W26
color orange, site_W26
label site_W26, "W26"
deselect
orient 1prg
# when the view is right:  ray 2400, 1800  then  png figures/sites_gamma.png, dpi=300
