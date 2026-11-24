# PPARα: conserved water sites (>= 30% at 1 Å) and nearby helices
# run from the project folder:  pymol figures/pymol_sites_alpha.pml
reinitialize
load data/raw/PPARalpha/1k7l.cif, 1k7l
remove 1k7l and not chain A
bg_color white
hide everything
show cartoon, 1k7l
color grey90, 1k7l
set cartoon_transparency, 0.15
set label_size, 16
set label_color, black
set float_labels, 1
select hx_H1, 1k7l and resi 202+203+204+205+206+207+208+209+210+211+212+213+214+215+216+217
color paleyellow, hx_H1
pseudoatom lab_H1, hx_H1 and name CA
label lab_H1, "H1"
hide everything, lab_H1
show labels, lab_H1
select hx_H2, 1k7l and resi 222+223+224+225+226+227+228+229
color grey70, hx_H2
pseudoatom lab_H2, hx_H2 and name CA
label lab_H2, "H2"
hide everything, lab_H2
show labels, lab_H2
select hx_H2p, 1k7l and resi 244+245+246+247+248+249+250+251
color grey60, hx_H2p
pseudoatom lab_H2p, hx_H2p and name CA
label lab_H2p, "H2'"
hide everything, lab_H2p
show labels, lab_H2p
select hx_H3, 1k7l and resi 267+268+269+270+271+272+273+274+275+276+277+278+279+280+281+282+283+284+285+286+287+288+289+290+291+292+293+294
color slate, hx_H3
pseudoatom lab_H3, hx_H3 and name CA
label lab_H3, "H3"
hide everything, lab_H3
show labels, lab_H3
select hx_H4_5, 1k7l and resi 297+298+299+300+301+302+303+304+305+306+307+308+309+310+311+312+313+314+315+316+317+318+319+320+321+322+323+324
color aquamarine, hx_H4_5
pseudoatom lab_H4_5, hx_H4_5 and name CA
label lab_H4_5, "H4/5"
hide everything, lab_H4_5
show labels, lab_H4_5
select hx_H6, 1k7l and resi 341+342+343+344+345+346
color wheat, hx_H6
pseudoatom lab_H6, hx_H6 and name CA
label lab_H6, "H6"
hide everything, lab_H6
show labels, lab_H6
select hx_H7, 1k7l and resi 350+351+352+353+354+355+356+357+358+359+360+361+362+363+364+365+366+367+368
color lightorange, hx_H7
pseudoatom lab_H7, hx_H7 and name CA
label lab_H7, "H7"
hide everything, lab_H7
show labels, lab_H7
select hx_H8, 1k7l and resi 372+373+374+375+376+377+378+379+380+381+382+383
color palegreen, hx_H8
pseudoatom lab_H8, hx_H8 and name CA
label lab_H8, "H8"
hide everything, lab_H8
show labels, lab_H8
select hx_H9, 1k7l and resi 394+395+396+397+398+399+400+401+402+403+404+405+406+407+408+409+410+411+412+413+414+415
color lightpink, hx_H9
pseudoatom lab_H9, hx_H9 and name CA
label lab_H9, "H9"
hide everything, lab_H9
show labels, lab_H9
select hx_H10_11, 1k7l and resi 422+423+424+425+426+427+428+429+430+431+432+433+434+435+436+437+438+439+440+441+442+443+444+445+446+447+448+449+450
color violet, hx_H10_11
pseudoatom lab_H10_11, hx_H10_11 and name CA
label lab_H10_11, "H10/11"
hide everything, lab_H10_11
show labels, lab_H10_11
select hx_H12, 1k7l and resi 457+458+459+460+461+462+463+464+465+466+468
color salmon, hx_H12
pseudoatom lab_H12, hx_H12 and name CA
label lab_H12, "H12"
hide everything, lab_H12
show labels, lab_H12
select site_W1, 1k7l and resn HOH and resi 488
show spheres, site_W1
set sphere_scale, 0.45, site_W1
color red, site_W1
label site_W1, "W1"
select site_W2, 1k7l and resn HOH and resi 494
show spheres, site_W2
set sphere_scale, 0.45, site_W2
color red, site_W2
label site_W2, "W2"
select site_W3, 1k7l and resn HOH and resi 492
show spheres, site_W3
set sphere_scale, 0.45, site_W3
color red, site_W3
label site_W3, "W3"
select site_W4, 1k7l and resn HOH and resi 515
show spheres, site_W4
set sphere_scale, 0.3, site_W4
color orange, site_W4
label site_W4, "W4"
select site_W5, 1k7l and resn HOH and resi 487
show spheres, site_W5
set sphere_scale, 0.3, site_W5
color orange, site_W5
label site_W5, "W5"
select site_W6, 1k7l and resn HOH and resi 482
show spheres, site_W6
set sphere_scale, 0.3, site_W6
color orange, site_W6
label site_W6, "W6"
select site_W7, 1k7l and resn HOH and resi 477
show spheres, site_W7
set sphere_scale, 0.3, site_W7
color orange, site_W7
label site_W7, "W7"
select site_W8, 1k7l and resn HOH and resi 507
show spheres, site_W8
set sphere_scale, 0.3, site_W8
color orange, site_W8
label site_W8, "W8"
select site_W9, 1k7l and resn HOH and resi 493
show spheres, site_W9
set sphere_scale, 0.3, site_W9
color orange, site_W9
label site_W9, "W9"
select site_W10, 1k7l and resn HOH and resi 491
show spheres, site_W10
set sphere_scale, 0.3, site_W10
color orange, site_W10
label site_W10, "W10"
select site_W11, 1k7l and resn HOH and resi 512
show spheres, site_W11
set sphere_scale, 0.3, site_W11
color orange, site_W11
label site_W11, "W11"
deselect
orient 1k7l
# when the view is right:  ray 2400, 1800  then  png figures/sites_alpha.png, dpi=300
