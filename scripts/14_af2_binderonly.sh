#!/usr/bin/env bash
# Phase IG-mono (central campaign) — AlphaFold2 on the BINDER ALONE (foldability /
# self-consistency filter). For every MPNN sequence: thread it onto its binder backbone,
# emit binder-only (chain A), batch-predict with af2_initial_guess in monomer mode, and
# score by plddt_binder (fold confidence) + binder_aligned_rmsd (does it fold to the design?).
# This does NOT test binding — that's the later complex step. It just prunes 1600 -> shortlist.
#
# Input : outputs/13_mpnn_central/seqs/design_*.fa  +  outputs/11_rfd_central/design_*.pdb
# Output: outputs/14_af2_binderonly/{in,out}/ , scores.sc , ranked_binderonly.csv
# Usage : [GPU=0] bash scripts/14_af2_binderonly.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/common.sh"

GPU="${GPU:-0}"
SEQDIR="$PROJECT/outputs/13_mpnn_central/seqs"
BBDIR="$PROJECT/outputs/11_rfd_central"
OUT="$PROJECT/outputs/14_af2_binderonly"
INDIR="$OUT/in"; PREDOUT="$OUT/out"; SC="$OUT/scores.sc"
mkdir -p "$INDIR" "$PREDOUT"

echo "→ Threading MPNN sequences onto binder backbones (binder-only)..."
"$IG_PY" - "$SEQDIR" "$BBDIR" "$INDIR" "$BINDER_CHAIN" <<'PY'
import sys, glob, os
from pathlib import Path
ONE_TO_THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
"G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
"S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
seqdir, bbdir, indir, bchain = sys.argv[1:5]
def read_seqs(fa):  # skip record 0 (native poly-G); return designed sequences in order
    recs=[]; cur=None
    for line in open(fa):
        if line.startswith(">"):
            cur=[]; recs.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    seqs=["".join(r) for r in recs]
    return seqs[1:] if len(seqs)>1 else []
n=0
for fa in sorted(glob.glob(os.path.join(seqdir,"design_*.fa"))):
    name=os.path.basename(fa)[:-3]
    bb=os.path.join(bbdir,name+".pdb")
    if not os.path.exists(bb): continue
    atoms=[l for l in open(bb) if l.startswith("ATOM") and l[21]==bchain]
    order=[]; seen=set()
    for l in atoms:
        if l[22:26] not in seen: seen.add(l[22:26]); order.append(l[22:26])
    ridx={r:i for i,r in enumerate(order)}
    for i,seq in enumerate(read_seqs(fa)):
        if len(seq)!=len(order): continue
        out=[]; ser=0
        for l in atoms:
            ri=ridx[l[22:26]]; ser+=1
            out.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[ri]]:>3s} A{ri+1:4d}{l[26:]}")
        out+=["TER","END"]
        Path(os.path.join(indir,f"{name}_s{i}.pdb")).write_text("\n".join(out)+"\n")
        n+=1
print(f"  threaded {n} binder-only structures")
PY

NIN=$(find "$INDIR" -name '*.pdb' | wc -l)
echo "→ Running af2_initial_guess (monomer mode) over $NIN structures (one model load)..."
cd "$IG_DIR"
CUDA_VISIBLE_DEVICES=$GPU "$IG_PY" predict.py \
  -pdbdir "$INDIR" -outpdbdir "$PREDOUT" -scorefilename "$SC" \
  -checkpoint_name "$OUT/check.point"

echo "→ Ranking by fold confidence (plddt) + self-consistency (rmsd)..."
"$IG_PY" - "$SC" "$PROJECT/outputs/11_rfd_central/monomer_contacts.csv" "$OUT/ranked_binderonly.csv" <<'PY'
import sys, csv, re
sc, geomcsv, outcsv = sys.argv[1:4]
# geometry: monomers-contacted per backbone
geom={}
try:
    for r in csv.DictReader(open(geomcsv)):
        geom[r["design"]]={"n_monomers":int(r["n_monomers"]),"weakest_chain":int(r["weakest_chain"]),
                           "binder_rxy":float(r["binder_rxy"]),"exterior_cap":r["exterior_cap"]}
except FileNotFoundError: pass
rows=[]
hdr=None
for line in open(sc):
    if not line.startswith("SCORE:"): continue
    t=line.split()[1:]
    if t[0]=="binder_aligned_rmsd": hdr=t; continue
    if hdr and len(t)==len(hdr):
        d=dict(zip(hdr,t)); tag=re.sub(r"_af2pred$","",d["description"])
        bb=re.match(r"(design_\d+)_s\d+",tag)
        bb=bb.group(1) if bb else tag
        g=geom.get(bb,{})
        rows.append({"tag":tag,"backbone":bb,
            "plddt_binder":round(float(d["plddt_binder"]),2),
            "binder_rmsd":round(float(d["binder_aligned_rmsd"]),2),
            "n_monomers":g.get("n_monomers",""),"weakest_chain":g.get("weakest_chain",""),
            "binder_rxy":g.get("binder_rxy",""),"exterior_cap":g.get("exterior_cap",""),
            "foldable":bool(float(d["plddt_binder"])>80 and float(d["binder_aligned_rmsd"])<2.0)})
# rank: foldable first, then plddt desc, then rmsd asc
rows.sort(key=lambda r:(not r["foldable"], -r["plddt_binder"], r["binder_rmsd"]))
keys=["rank","tag","backbone","plddt_binder","binder_rmsd","foldable","n_monomers","weakest_chain","binder_rxy","exterior_cap"]
with open(outcsv,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=keys,lineterminator="\n"); w.writeheader()
    for i,r in enumerate(rows,1): w.writerow({"rank":i,**r})
nf=sum(r["foldable"] for r in rows)
print(f"  {nf}/{len(rows)} foldable (plddt>80 AND rmsd<2A) -> {outcsv}")
for r in rows[:8]:
    print(f"  {r['tag']}: plddt={r['plddt_binder']} rmsd={r['binder_rmsd']} all6={r['n_monomers']==6} foldable={r['foldable']}")
PY
echo "→ Done. Ranked: $OUT/ranked_binderonly.csv"
