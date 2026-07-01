#!/usr/bin/env python3
"""Refined T=1 vs T=3: partner-pairing. A chain-A residue r is T=3-UNIQUE if it contacts a partner
residue s in the T=3 hexamer that it NEVER contacts anywhere in the full T=1 capsid (1LP3, 60-mer).
This catches seams that are 'at an interface in both' but with different partners/geometry.
Updates outputs/t1t3_site_table.csv (adds t3_unique_partner) and regenerates the combined figure.
"""
from __future__ import annotations
import csv, numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner

ROOT="/data/binder_software/pre-binder"
HEX=f"{ROOT}/inputs/151lp3t3_hexamer_6chain.pdb"; T1=f"{ROOT}/inputs/1LP3_T1_aav2.pdb"
RMSF=f"{ROOT}/inputs/rmsf_BB_COM_atomcenter.xvg"; CUT=5.0
T2O={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
 'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def car(model,cid):
    xyz,rid=[],[]
    for res in model[cid]:
        if res.resname not in T2O: continue
        for a in res:
            if a.element!="H": xyz.append(a.coord); rid.append(res.id[1])
    return np.array(xyz), np.array(rid)
def seq(model,cid):
    it=[(r.id[1],T2O[r.resname]) for r in model[cid] if r.resname in T2O]
    return "".join(a for _,a in it),[n for n,_ in it]
def biomt(fn):
    ops={}
    for l in open(fn):
        if l[11:23].strip().startswith("BIOMT"):
            s=l.split(); ops.setdefault(int(s[3]),np.zeros((3,4)))[int(s[2][-1])-1]=list(map(float,s[4:8]))
    return [ops[k] for k in sorted(ops)]

p=PDBParser(QUIET=True); hexm=p.get_structure("h",HEX)[0]; t1m=p.get_structure("t",T1)[0]
# numbering map 1LP3 -> hex
sH,nH=seq(hexm,"A"); sT,nT=seq(t1m,"A")
al=PairwiseAligner(); al.mode="global"; al.open_gap_score=-5; al.extend_gap_score=-0.5
a=al.align(sH,sT)[0]; m_t1_hex={}
for (h0,h1),(t0,t1b) in zip(a.aligned[0],a.aligned[1]):
    for k in range(h1-h0): m_t1_hex[nT[t0+k]]=nH[h0+k]

# T=3 partner residues (hex numbering)
aX,aR=car(hexm,"A"); oX=[];oR=[]
for c in "BCDEF":
    x,r=car(hexm,c); oX.append(x); oR.append(r)
oX=np.vstack(oX); oR=np.concatenate(oR); ot=cKDTree(oX)
t3p={}
for x,r in zip(aX,aR):
    for j in ot.query_ball_point(x,CUT): t3p.setdefault(int(r),set()).add(int(oR[j]))

# T=1 partner residues over all 60 images (1LP3 numbering) -> hex
tX,tR=car(t1m,"A"); rt=cKDTree(tX); t1p={}
for op in biomt(T1):
    img=(tX@op[:,:3].T)+op[:,3]
    if np.allclose(img,tX,atol=1e-3): continue
    pairs=rt.query_ball_tree(cKDTree(img),CUT)
    for i,h in enumerate(pairs):
        if h:
            r=int(tR[i])
            for j in h: t1p.setdefault(r,set()).add(int(tR[j]))
t1p_hex={}
for r1,ps in t1p.items():
    if r1 in m_t1_hex:
        t1p_hex.setdefault(m_t1_hex[r1],set()).update(m_t1_hex[s] for s in ps if s in m_t1_hex)

# RMSF avg per residue (nm->A)
racc={}
for l in open(RMSF):
    l=l.strip()
    if l and l[0] not in "#@":
        r,v=l.split()[:2]; racc.setdefault(int(r),[]).append(float(v)*10)
rmsf={r:float(np.mean(v)) for r,v in racc.items()}

site={int(r["resid"]):r for r in csv.DictReader(open(f"{ROOT}/outputs/site_analysis.csv"))}
rows=[]
for rid in sorted(site):
    s=site[rid]; T3=t3p.get(rid,set()); T1=t1p_hex.get(rid,set())
    uniq=T3-T1
    rows.append({"resid":rid,"aa":s["aa"],"rsasa":float(s["rsasa"]),"phi":float(s["phi"]),
        "concavity":int(s["concavity"]),"n_chains":int(s["n_chains"]),"rmsf":round(rmsf.get(rid,np.nan),2),
        "t3_partners":len(T3),"t1_partners":len(T1),"t3_unique_partners":len(uniq),
        "t3_unique":bool(len(uniq)>0)})
dthr=np.percentile([r["concavity"] for r in rows if r["rsasa"]>0.25],60)
pthr=np.percentile([abs(r["phi"]) for r in rows],60)
for r in rows:
    r["candidate"]=bool(r["rsasa"]>0.25 and r["t3_unique"] and r["concavity"]>=dthr and abs(r["phi"])>=pthr)
with open(f"{ROOT}/outputs/t1t3_site_table.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),lineterminator="\n"); w.writeheader(); w.writerows(rows)

uniq=sum(r["t3_unique"] for r in rows); eu=sum(r["t3_unique"] and r["rsasa"]>0.25 for r in rows)
cand=[r["resid"] for r in rows if r["candidate"]]
sites,cur=[],[]
for rid in cand:
    if cur and rid-cur[-1]>3: sites.append(cur); cur=[]
    cur.append(rid)
if cur: sites.append(cur)
print(f"PARTNER-PAIRING: T=3-unique residues={uniq} (exposed={eu}); candidates={len(cand)} in {len(sites)} seams")
for st in sorted(sites,key=len,reverse=True)[:10]:
    m=[r for r in rows if r["resid"] in st]; ph=np.mean([r["phi"] for r in m])
    print(f"  {st[0]}-{st[-1]} ({len(st)} res) phi={ph:+.2f}({'acid' if ph<0 else 'base'}) "
          f"rSASA={np.mean([r['rsasa'] for r in m]):.2f} rmsf={np.nanmean([r['rmsf'] for r in m]):.1f}")

# B-factor PDB for ChimeraX (T=3-unique-partner count)
um={r["resid"]:r["t3_unique_partners"] for r in rows}
out=[]
for l in open(HEX):
    if l.startswith("ATOM") and l[21]=="A":
        out.append(f"{l[:60]}{float(um.get(int(l[22:26]),0)):6.2f}{l[66:].rstrip()}")
open(f"{ROOT}/outputs/chainA_t3unique.pdb","w").write("\n".join(out)+"\n")

# combined figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Liberation Sans","Arial","DejaVu Sans"],
 "font.size":8,"axes.titlesize":9.5,"axes.labelsize":8.5,"xtick.labelsize":7,"ytick.labelsize":7,
 "legend.fontsize":7,"text.color":"black","axes.labelcolor":"black","axes.titlecolor":"black",
 "xtick.color":"black","ytick.color":"black","axes.edgecolor":"black","pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.8})
x=np.array([r["resid"] for r in rows]); G=lambda k:np.array([r[k] for r in rows],float)
fig,ax=plt.subplots(6,1,figsize=(7.2,9.4),sharex=True)
fig.suptitle("Site selection for a T=3-specific binder (partner-pairing T=1 vs T=3)",fontsize=10.5,fontweight="bold",y=0.998)
def deco(a):
    for st in sites: a.axvspan(st[0],st[-1],color="#fff2a8",alpha=0.95,lw=0)
ax[0].plot(x,G("rsasa"),color="#0b7a6f",lw=.9); ax[0].set_ylabel("Rel. SASA"); ax[0].set_ylim(0,1)
ax[1].fill_between(x,0,G("phi"),where=(G("phi")>=0),color="#3b6ea5",alpha=.85,lw=0)
ax[1].fill_between(x,0,G("phi"),where=(G("phi")<0),color="#c0504d",alpha=.85,lw=0); ax[1].axhline(0,color="k",lw=.5)
ax[1].set_ylabel("Electrostatic φ")
ax[2].plot(x,G("concavity"),color="#c4622d",lw=.9); ax[2].set_ylabel("Concavity")
ax[3].plot(x,G("rmsf"),color="#7b5ea7",lw=.9); ax[3].set_ylabel("MD RMSF (Å)")
ax[4].plot(x,G("t3_partners"),color="#1d3557",lw=.9,label="T=3"); ax[4].plot(x,G("t1_partners"),color="#e07b39",lw=.9,label="T=1")
ax[4].set_ylabel("# partner res"); ax[4].legend(frameon=False,loc="upper right",ncol=2)
ax[5].fill_between(x,0,G("t3_unique_partners"),color="#2a7a2a",lw=0,step="mid"); ax[5].set_ylabel("T=3-UNIQUE\npartners")
for a in ax: deco(a); a.margins(x=0); a.tick_params(width=0.8,length=2.5)
ax[-1].set_xlabel("VP3 residue number")
fig.text(0.5,0.004,"Yellow = recommended spots (exposed+concave+charged+T=3-unique). Panel 6: number of partner residues "
 "contacted in T=3 that are NEVER contacted in the full T=1 capsid (genuine assembly-unique contacts).",ha="center",fontsize=7)
fig.tight_layout(rect=[0,0.012,1,0.985])
for e in ("pdf","png","svg"): fig.savefig(f"{ROOT}/reports/site_selection_combined.{e}",dpi=300)
print("wrote reports/site_selection_combined.{pdf,png,svg}; outputs/t1t3_site_table.csv; outputs/chainA_t3unique.pdb")
