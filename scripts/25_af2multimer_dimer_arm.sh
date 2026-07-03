#!/usr/bin/env bash
# Stage 1 Step 5 — AF2 (initial-guess) COMPLEX cross-check of the Rosetta hits.
# Orthogonal to Rosetta: does AF2 predict a confident binder<->dimer interface? Rejects the
# greasy-collapse Rosetta artifacts (very negative dG but high buns). For each hit: build a 2-chain
# complex (binder -> chain A threaded with its MPNN seq; target A+F -> merged chain B), then run
# af2_initial_guess and read pae_interaction (want <10) + plddt_binder (want >80).
#
# Input : tag list (arg 1, default = the 9 neg-dG + Lys42-saltbridge hits) + 01_rfd/design_*.pdb
# Output: <campaign>/05_af2multimer/{in,out}/ , scores.sc , ranked_complex.csv
# Usage : [GPU=0] bash scripts/25_af2multimer_dimer_arm.sh [tags.txt]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; source "$HERE/lib/common.sh"
GPU="${GPU:-0}"; BCHAIN="${BCHAIN:-B}"
TARGET_SHELL="${TARGET_SHELL:-0}"   # 0 = full A+F target; >0 = keep only target residues within N A of binder
CAMP="$PROJECT/outputs/C3_symmetric_Binder_2026_07_02"
RFDDIR="$CAMP/01_rfd"; SEQDIR="$CAMP/02_mpnn/seqs"
if [ "$TARGET_SHELL" = "0" ]; then OUT="$CAMP/05_af2multimer"; else OUT="$CAMP/05_af2multimer/trim${TARGET_SHELL}"; fi
INDIR="$OUT/in"; PREDOUT="$OUT/out"; SC="$OUT/scores.sc"
mkdir -p "$INDIR" "$PREDOUT"
TAGS="${1:-$CAMP/04_rosetta/af2mm_tags.txt}"
echo "→ TARGET_SHELL=$TARGET_SHELL  (0=full target; >0=trim target to that shell around binder)"

echo "→ Building 2-chain complexes (binder→A threaded, target A+F→merged B) for tags in $TAGS ..."
"$IG_PY" - "$TAGS" "$RFDDIR" "$SEQDIR" "$INDIR" "$BCHAIN" "$TARGET_SHELL" <<'PY'
import sys, os
from pathlib import Path
ONE_TO_THREE={"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
"G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO",
"S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
tags_f, rfddir, seqdir, indir, bchain, shell = sys.argv[1:7]
shell=float(shell)
def seq_for(tag):
    bb,s=tag.rsplit("_s",1); recs,cur=[],None
    for line in open(f"{seqdir}/{bb}.fa"):
        if line.startswith(">"): cur=[]; recs.append(cur)
        elif cur is not None: cur.append(line.strip())
    return ["".join(r) for r in recs][1:][int(s)]
def xyz(l): return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
def d2(a,b): return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
n=0
for tag in [t.strip() for t in open(tags_f) if t.strip()]:
    bb=tag.rsplit("_s",1)[0]; design=f"{rfddir}/{bb}.pdb"
    lines=open(design).read().splitlines()
    binder=[l for l in lines if l.startswith("ATOM") and l[21]==bchain]
    bxyz=[xyz(l) for l in binder]
    order,seen=[],set()
    for l in binder:
        if l[22:26] not in seen: seen.add(l[22:26]); order.append(l[22:26])
    seq=seq_for(tag); ridx={r:i for i,r in enumerate(order)}
    if len(seq)!=len(order):
        print(f"  skip {tag}: seq {len(seq)} != binder {len(order)}"); continue
    out=[]; ser=0
    # binder -> chain A, threaded
    for l in binder:
        i=ridx[l[22:26]]; ser+=1
        out.append(f"{l[:6]}{ser:5d}{l[11:17]}{ONE_TO_THREE[seq[i]]:>3s} A{i+1:4d}{l[26:]}")
    out.append("TER")
    # target chains A + F -> merged chain B, continuous numbering (physical gaps => AF2 chainbreaks).
    # If shell>0, keep only target residues with any atom within `shell` A of the binder (so AF2 can
    # hold a small rigid epitope instead of drifting the full 1008-res target).
    rnum=0; sh2=shell*shell; kept=0
    for c in ("A","F"):
        # group this chain's atoms by residue
        resatoms={}
        for l in lines:
            if l.startswith("ATOM") and l[21]==c:
                resatoms.setdefault(l[22:26],[]).append(l)
        for rid,atoms in resatoms.items():
            if shell>0:
                near=any(d2(xyz(a),bx)<sh2 for a in atoms for bx in bxyz)
                if not near: continue
            rnum+=1; kept+=1
            for l in atoms:
                ser+=1
                out.append(f"{l[:6]}{ser:5d}{l[11:21]}B{rnum:4d}{l[26:]}")
    out+=["TER","END"]
    Path(f"{indir}/{tag}.pdb").write_text("\n".join(out)+"\n"); n+=1
    if shell>0: print(f"  {tag}: kept {kept} target residues within {shell:.0f} A of binder")
print(f"  built {n} complexes")
PY

NIN=$(find "$INDIR" -name '*.pdb' | wc -l)
echo "→ af2_initial_guess (complex mode) over $NIN structures on GPU$GPU ..."
cd "$IG_DIR"
CUDA_VISIBLE_DEVICES=$GPU "$IG_PY" predict.py \
  -pdbdir "$INDIR" -outpdbdir "$PREDOUT" -scorefilename "$SC" -checkpoint_name "$OUT/check.point"

echo "→ Ranking by pae_interaction (pass = pae_interaction<10 AND plddt_binder>80)..."
"$IG_PY" - "$SC" "$CAMP/04_rosetta/relaxed_metrics.csv" "$OUT/ranked_complex.csv" <<'PY'
import sys, csv, re
sc, relcsv, outcsv = sys.argv[1:4]
rel={r["tag"]:r for r in csv.DictReader(open(relcsv))}
rows=[]; hdr=None
for line in open(sc):
    if not line.startswith("SCORE:"): continue
    t=line.split()[1:]
    if not hdr and "pae_interaction" in t: hdr=t; continue
    if hdr and len(t)==len(hdr):
        d=dict(zip(hdr,t)); tag=re.sub(r"_af2pred$","",d["description"])
        rr=rel.get(tag,{})
        rows.append({"tag":tag,
            "pae_interaction":round(float(d["pae_interaction"]),2),
            "plddt_binder":round(float(d["plddt_binder"]),2),
            "target_rmsd":round(float(d.get("target_aligned_rmsd",d.get("binder_aligned_rmsd",0))),2),
            "dG_relaxed":rr.get("dG_relaxed",""),"buns":rr.get("delta_unsat_hbonds",""),
            "K42sb":rr.get("salt_bridge_K42",""),
            "pass":bool(float(d["pae_interaction"])<10 and float(d["plddt_binder"])>80)})
rows.sort(key=lambda r:r["pae_interaction"])
keys=["tag","pae_interaction","plddt_binder","target_rmsd","dG_relaxed","buns","K42sb","pass"]
with open(outcsv,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=keys,extrasaction="ignore",lineterminator="\n")
    w.writeheader(); w.writerows(rows)
np=sum(r["pass"] for r in rows)
print(f"  {np}/{len(rows)} PASS (pae_interaction<10 AND plddt_binder>80) -> {outcsv}")
for r in rows:
    print(f"  {r['tag']}: pae_i={r['pae_interaction']} plddt={r['plddt_binder']} "
          f"dG={r['dG_relaxed']} buns={r['buns']} K42sb={r['K42sb']} pass={r['pass']}")
PY
echo "→ Done. Ranked: $OUT/ranked_complex.csv"
