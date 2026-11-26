"""Helix assignment and contact distances of the top reference waters (1K7L, 1PRG, 1GWX).
Secondary-structure boundaries: supervisor's table (PPAR LBD Secondary Structure Boundaries)."""
import sys, numpy as np
from Bio.PDB import MMCIFParser
SS = {  # element: (gamma1, alpha, delta)
 "H1":((208,226),(199,217),(172,190)), "H2":((230,238),(221,229),(194,202)),
 "B1":((246,250),(237,241),(210,214)), "H2p":((252,261),(243,252),(216,225)),
 "H3":((274,300),(265,291),(238,264)), "H3p":((306,309),(297,300),(270,273)),
 "H4":((311,320),(302,311),(275,284)), "H5":((322,333),(313,324),(286,297)),
 "B2":((334,338),(325,329),(298,302)), "B3":((339,342),(330,333),(303,306)),
 "B4":((346,349),(337,340),(310,313)), "H6":((350,355),(341,346),(314,319)),
 "H7":((361,378),(352,369),(325,342)), "H8":((380,394),(371,385),(344,358)),
 "H9":((403,424),(394,415),(367,388)), "H10":((430,457),(421,448),(394,421)),
 "H12":((467,475),(458,466),(431,439))}
SETS = {"gamma":("1prg","data/raw/PPARgamma/1prg.cif",0,{"W1":489,"W2":614,"W3":537,"W4":498,"W5":479}),
        "alpha":("1k7l","data/raw/PPARalpha/1k7l.cif",1,{"W1":488,"W2":494,"W3":492}),
        "delta":("1gwx","data/raw/PPARdelta/1gwx.cif",2,{"W1":514,"W2":511,"W3":485})}
OFFSET = {0:0, 1:0, 2:36}   # 1GWX is numbered +36 vs the table (H3 Cys249 = Cys285 in 1GWX)
def element(resi, col):
    for k,v in SS.items():
        a,b = v[col][0]+OFFSET[col], v[col][1]+OFFSET[col]
        if a <= resi <= b: return k
    return "loop"
base = sys.argv[1]
rows=[]
for sub,(pdb,path,col,ws) in SETS.items():
    s = MMCIFParser(QUIET=True).get_structure(pdb, f"{base}/{path}")[0]["A"]
    prot = [a for r in s if r.id[0]==" " for a in r]
    waters = {r.id[1]: r for r in s if r.id[0]=="W"}
    otherw = [a for r in waters.values() for a in r]
    print(f"\n== PPAR{sub} {pdb} chain A: residues {min(r.id[1] for r in s if r.id[0]==' ')}-{max(r.id[1] for r in s if r.id[0]==' ')}")
    for w,resi in ws.items():
        O = waters[resi]["O"]
        d = [(np.linalg.norm(a.coord-O.coord),a) for a in prot]
        d.sort(key=lambda x:x[0])
        near_el = {}
        for dist,a in d:
            if dist>4.0: break
            el = element(a.get_parent().id[1],col); near_el[el]=min(near_el.get(el,9),dist)
        pol = [(dist,a) for dist,a in d if dist<=3.5 and a.element in ("N","O")]
        wat = sorted((np.linalg.norm(a.coord-O.coord),a.get_parent().id[1]) for a in otherw if a.get_parent().id[1]!=resi)
        wat = [x for x in wat if x[0]<=3.5]
        print(f"{w} HOH {resi}  B={O.bfactor:.1f}  nearest elements(<=4A): " + ", ".join(f"{k} {v:.2f}" for k,v in sorted(near_el.items(), key=lambda x:x[1])))
        for dist,a in pol:
            r=a.get_parent()
            print(f"    {r.get_resname()}{r.id[1]}:{a.get_id()} [{element(r.id[1],col)}] {dist:.2f}")
            rows.append((sub,pdb,w,resi,r.get_resname(),r.id[1],a.get_id(),element(r.id[1],col),round(float(dist),2)))
        if not pol: print(f"    no protein N/O within 3.5 A; nearest atom {d[0][1].get_parent().get_resname()}{d[0][1].get_parent().id[1]}:{d[0][1].get_id()} {d[0][0]:.2f}")
        for dist,r in wat: print(f"    water HOH{r} {dist:.2f}")
import csv
with open(f"{base}/data/results/top_water_contacts.csv","w",newline="") as f:
    wr=csv.writer(f); wr.writerow(["subtype","pdb","site","water_resi","residue","resi","atom","element","distance_A"]); wr.writerows(rows)
